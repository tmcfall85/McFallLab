import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM
import sys


def main(prompt_in):
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print("CUDA is available. Using GPU.")
    else:
        device = torch.device("cpu")
        print("CUDA is not available. Using CPU.")

    tokenizer = AutoTokenizer.from_pretrained("vandijklab/C2S-Scale-Gemma-2-2B")
    model = AutoModelForCausalLM.from_pretrained("vandijklab/C2S-Scale-Gemma-2-2B")
    model.to(device)
    prompt_in_filename = f"{prompt_in}.csv"
    prompts = pd.read_csv(prompt_in_filename)
    out_filename = prompt_in_filename.replace("prompts", "embeddings")

    all_embeddings = []
    all_acc_ids = []
    all_layer_ids = []
    cols = []
    for i in range(len(prompts)):
        prompt = prompts.iloc[i].prompt
        accession_id = prompts.iloc[i].accession_id
        print(f"Accession id: {accession_id}")
        inputs = tokenizer(prompt, return_tensors="pt")
        inputs.to(device)
        # next time do max_new_tokens = 2 nd remove range(10) part and just do [-1][-1][-1]
        generate_ids = model.generate(
            inputs.input_ids,
            max_new_tokens=1,
            return_dict_in_generate=True,
            output_hidden_states=True,
        )

        stacked_tensors = torch.stack(generate_ids["hidden_states"][-1])

        for i in range(20):
            all_embeddings.append(
                stacked_tensors[-1][-1][(i + 1) * -1].cpu().detach().numpy()
            )
            all_acc_ids.append(accession_id)
            all_layer_ids.append((i + 1) * -1)

        df_so_far = pd.DataFrame(all_embeddings, index=[all_acc_ids, all_layer_ids])
        df_so_far.index.set_names(["accession_id", "layer_id"], inplace=True)
        df_so_far.to_csv(out_filename)
        print(f"saving df snapshot, len{len(df_so_far)}")
    print(f"done: {out_filename.split('.csv')[0]}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_modeling.py <prompt_in>")
        sys.exit(1)
    main(sys.argv[1])
