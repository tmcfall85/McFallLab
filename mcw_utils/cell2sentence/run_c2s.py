import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM
import sys
from pathlib import Path
from datetime import datetime
import json


def main(prompt_folder):

    results_dir = Path(f"embeddings/{str(datetime.now()).replace(' ','_')}")
    results_dir.mkdir(exist_ok=True)
    for file_name in prompt_folder.iterdir():
        if file_name.is_file():
            if file_name.suffix == ".json":
                print(file_name)
                with open(file_name, "r") as fp:
                    metadata = json.load(fp)
                train_embeddings_filenames = []
                holdout_embeddings_filenames = []
                
                for prompt_train_fname in metadata["prompt_train_fnames"]:
                    print(f"Running c2s on {prompt_train_fname}")
                    train_embeddings_filename = c2s(prompt_train_fname, results_dir)
                    train_embeddings_filenames.append(str(train_embeddings_filename))

                metadata["train_embedding_filenames"] = train_embeddings_filenames
                
                for prompt_holdout_fname in metadata["prompt_holdout_fnames"]:
                    print(f"Running c2s on {prompt_holdout_fname}")
                    holdout_embeddings_filename = c2s(prompt_holdout_fname, results_dir)
                    holdout_embeddings_filenames.append(str(holdout_embeddings_filename))

                metadata["holdout_embedding_filenames"] = holdout_embeddings_filenames

                with open(results_dir / f"embedding_{file_name.stem}.json", "w") as fp:
                    json.dump(metadata, fp)


def c2s(prompt_in, results_dir):
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print("CUDA is available. Using GPU.")
    else:
        device = torch.device("cpu")
        print("CUDA is not available. Using CPU.")

    tokenizer = AutoTokenizer.from_pretrained("vandijklab/C2S-Scale-Gemma-2-2B")
    model = AutoModelForCausalLM.from_pretrained("vandijklab/C2S-Scale-Gemma-2-2B")
    model.to(device)
    # prompt_in_filename = f"{prompt_in}.csv"
    prompts = pd.read_csv(prompt_in)
    in_fname = Path(prompt_in)
    out_fname_stem = in_fname.stem
    out_filename = results_dir / f"embeddings_{out_fname_stem}.csv"
    # prompt_df_with_outcomes.to_csv(
    #    out_fname,
    #    index=False,
    # )
    # out_filename = prompt_in_filename.replace("prompts", "embeddings")

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
    print(f"done: {out_filename}")
    return out_filename


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_modeling.py <prompt_in>")
        sys.exit(1)
    main(Path(sys.argv[1]))
