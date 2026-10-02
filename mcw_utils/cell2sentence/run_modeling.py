import sys
from sklearn.manifold import TSNE
import plotly.express as px
import pandas as pd
from datetime import datetime
from pathlib import Path
import json
import pickle 

import numpy as np
from FCNet_prekfold import train_fc_network_prekfold
from FCNet import train_fc_network
from FCNet_single import train_fc_network_single
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt


def embedding_to_X_y(embeddings_in, prepare_in, results_dir, max_layer):

    embeddings_df = pd.read_csv(embeddings_in)
    embeddings_df.set_index(["accession_id", "layer_id"], inplace=True)
    label_data = pd.read_csv(prepare_in)

    early_mid_late_panc_data = label_data[
        label_data.recurrence_time_sur.isin(["early", "mid", "late"])
    ]

    tsne = TSNE(n_components=2, random_state=0, perplexity=5)

    proj = tsne.fit_transform(embeddings_df)
    fig = px.scatter(
        proj,
        x=0,
        y=1,
        color=embeddings_df.reset_index()["layer_id"],
        labels={"color": "layer_id"},
    )
    fig.write_html(results_dir / f"tsne_plot_{Path(embeddings_in).stem}.html")

    X = []
    y = []
    acc_ids = []
    for i in range(len(early_mid_late_panc_data.accession_id.values)):
        row = early_mid_late_panc_data.iloc[i]
        acc_ids.append(row.accession_id)
        embeddings = []
        for i in range(max_layer):
            embeddings.append(
                embeddings_df.loc[(row.accession_id, (i + 1) * -1)].values
            )
        if row.recurrence_time_sur == "early":
            y.append(0)
            X.append(embeddings)
        elif row.recurrence_time_sur == "mid":
            y.append(1)
            X.append(embeddings)
        elif row.recurrence_time_sur == "late":
            y.append(2)
            X.append(embeddings)
    X = np.array(X).astype(np.float32)
    y = np.array(y)
    return X, y, acc_ids


def model(top_n_sig, embeddings_in, prepare_train, results_dir, kfolds_in=None, holdout_embeddings_in=None, prepare_holdout=None):
    max_layer = 9
    if kfolds_in is not None:
        Xs_dict = {}
        ys_dict = {}
        Xs = []
        ys = []
        Xs_h = []
        ys_h = []
        acc_ids = []
        for kfold in embeddings_in.keys():
            kfold_results_dir = results_dir / f'kfold_{kfold}'
            
            kfold_results_dir.mkdir(exist_ok=True)
            X, y, acc_id = embedding_to_X_y(embeddings_in[kfold], prepare_train[kfold], kfold_results_dir, max_layer)
            Xs.append(X)
            ys.append(y)
            if holdout_embeddings_in is not None and prepare_holdout is not None:
                X_h, y_h, acc_id_h = embedding_to_X_y(holdout_embeddings_in[kfold], prepare_holdout[kfold], kfold_results_dir, max_layer)
                
                Xs_h.append(X_h)
                ys_h.append(y_h)
                
            acc_ids.append(acc_id)
        models, histories, all_preds, all_preds_proba, all_y_vals, all_acc_ids, all_holdout_preds, all_holdout_preds_proba, all_holdout_y_vals, all_holdout_acc_ids = (
            train_fc_network_prekfold(
                Xs,
                ys,
                acc_ids,
                kfolds_in,
                Xs_np_h = Xs_h,
                ys_np_h = ys_h,
                acc_id_h = acc_id_h,
                hidden_sizes=(512, 128),
                n_classes=3,
                epochs=15,
                batch_size=2,
                lr=1e-3,
            )
        )
    else:
        X, y, acc_ids = embedding_to_X_y(embeddings_in, prepare_train, results_dir, max_layer)
        models, histories, all_preds, all_preds_proba, all_y_vals, all_acc_ids = (
            train_fc_network(
                X,
                y,
                kfolds_in,
                acc_ids,
                hidden_sizes=(512, 128),
                n_classes=3,
                epochs=15,
                batch_size=2,
                lr=1e-3,
            )
        )

    val = []
    for history in histories:
        val.append(history["val_acc"][-1])
    print(f"Max layers: {max_layer}")
    print(f"Average final validation acc: {np.mean(val):.3g} +- {np.std(val):.3g}")

    for i in range(len(all_preds)):
        print(f"Confusion matrix for fold {i+1}:")
        cm = confusion_matrix(all_y_vals[i], all_preds[i])
        print(cm)

    combined_all_y_vals = []
    combined_all_preds = []
    for all_y_val in all_y_vals:
        for y_val in all_y_val:
            combined_all_y_vals.append(y_val)
    for all_preds in all_preds:
        for preds in all_preds:
            combined_all_preds.append(preds)

    print(f"Combined confusion matrix :")
    cm = confusion_matrix(combined_all_y_vals, combined_all_preds)
    print(cm)

    early_v_mid_late_actual = []
    for all_y_val in all_y_vals:
        for y_val in all_y_val:
            if y_val == 0:
                early_v_mid_late_actual.append(0)
            else:
                early_v_mid_late_actual.append(1)
    early_v_mid_late_pred = []
    for all_pred in all_preds_proba:
        for y_pred in all_pred:
            early_v_mid_late_pred.append(1 - y_pred[0])
    early_v_mid_late_acc_id = []
    for val_acc_ids in all_acc_ids:
        for val_acc_id in val_acc_ids:
            early_v_mid_late_acc_id.append(val_acc_id)
    df_pred_out = pd.DataFrame(
        {
            "accession_id": early_v_mid_late_acc_id,
            "pred": early_v_mid_late_pred,
            "actual": early_v_mid_late_actual,
        }
    )
    df_pred_out.to_csv(results_dir / f"predictions_with_acc_ids_top_n_sig_{top_n_sig}.csv", index=False)
    fpr, tpr, thresholds = roc_curve(early_v_mid_late_actual, early_v_mid_late_pred)
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color="blue", label="ROC curve")
    plt.plot([0, 1], [0, 1], color="red", linestyle="--", label="Random guess")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate (Recall)")
    plt.title("ROC Curve")
    plt.legend()
    plt.savefig(results_dir / f"roc_curve_top_n_sig_{top_n_sig}.png")
    ras = roc_auc_score(early_v_mid_late_actual, early_v_mid_late_pred)

    print(f"ROC AUC score ({embeddings_in}): {ras}")
    if holdout_embeddings_in is not None and prepare_holdout is not None:
        '''
        holdout_embeddings_df = pd.read_csv(holdout_embeddings_in)
        holdout_embeddings_df.set_index(["accession_id", "layer_id"], inplace=True)
        label_data_holdout = pd.read_csv(prepare_holdout)

        early_mid_late_panc_data_holdout = label_data_holdout[
            label_data_holdout.recurrence_time_sur.isin(["early", "mid", "late"])
        ]
        #early_mid_late_panc_data_holdout = early_mid_late_data_holdout[
        #    early_mid_late_data_holdout.is_panc == True
        #]
        X_holdout = []
        y_holdout = []
        for i in range(len(early_mid_late_panc_data_holdout.accession_id.values)):
            row = early_mid_late_panc_data_holdout.iloc[i]
            embeddings_holdout = []
            for i in range(max_layer):
                embeddings_holdout.append(
                    holdout_embeddings_df.loc[(row.accession_id, (i + 1) * -1)].values
                )
            if row.recurrence_time_sur == "early":
                y_holdout.append(0)
                X_holdout.append(embeddings_holdout)
            elif row.recurrence_time_sur == "mid":
                y_holdout.append(1)
                X_holdout.append(embeddings_holdout)
            elif row.recurrence_time_sur == "late":
                y_holdout.append(2)
                X_holdout.append(embeddings_holdout)
        X_holdout = np.array(X_holdout).astype(np.float32)
        y_holdout = np.array(y_holdout)

        model, history, preds, preds_proba, y_vals = train_fc_network_single(
            X,
            y,
            X_holdout,
            y_holdout,
            hidden_sizes=(512, 128),
            n_classes=3,
            epochs=15,
            batch_size=2,
            lr=1e-3,
        )
        '''
        holdout = []
        for history in histories:
            holdout.append(history["holdout_acc"][-1])
        print(f"Max layers: {max_layer}")
        print(f"Average final holdout acc: {np.mean(holdout):.3g} +- {np.std(holdout):.3g}")

        for i in range(len(all_holdout_preds)):
            print(f"Confusion matrix for fold {i+1}:")
            cm = confusion_matrix(all_holdout_y_vals[i], all_holdout_preds[i])
            print(cm)

        combined_all_holdout_y_vals = []
        combined_all_holdout_preds = []
        early_v_mid_late_acc_id = []
        for val_acc_ids in all_holdout_acc_ids:
            for val_acc_id in val_acc_ids:
                early_v_mid_late_acc_id.append(val_acc_id)
        for all_holdout_y_val in all_holdout_y_vals:
            for y_val in all_holdout_y_val:
                combined_all_holdout_y_vals.append(y_val)
        for all_holdout_preds in all_holdout_preds:
            for preds in all_holdout_preds:
                combined_all_holdout_preds.append(preds)
        df_combo_pred_out = pd.DataFrame(
                    {
                        "accession_id": early_v_mid_late_acc_id,
                        "combo_pred": combined_all_holdout_preds,
                        "actual": combined_all_holdout_y_vals,
                    }
                )
        
        df_pred_out_mean = df_combo_pred_out.groupby('accession_id').mean().reset_index()
        def fix(row):
            return np.round(row.combo_pred)
        df_pred_out_mean['voted'] = df_pred_out_mean.apply(fix, axis=1)
        df_pred_out_mean.to_csv(results_dir / f"holdout_predictions_voted_agg_with_acc_ids_top_n_sig_{top_n_sig}.csv", index=False)
        print(f"Combined confusion matrix :")
        cm = confusion_matrix(df_pred_out_mean.actual, df_pred_out_mean.voted)
        print(cm)

        early_v_mid_late_actual = []
        for all_holdout_y_val in all_holdout_y_vals:
            for y_val in all_holdout_y_val:
                if y_val == 0:
                    early_v_mid_late_actual.append(0)
                else:
                    early_v_mid_late_actual.append(1)
        early_v_mid_late_pred = []
        for all_holdout_pred in all_holdout_preds_proba:
            for y_pred in all_holdout_pred:
                early_v_mid_late_pred.append(1 - y_pred[0])
        early_v_mid_late_acc_id = []
        for val_acc_ids in all_holdout_acc_ids:
            for val_acc_id in val_acc_ids:
                early_v_mid_late_acc_id.append(val_acc_id)
        df_pred_out = pd.DataFrame(
            {
                "accession_id": early_v_mid_late_acc_id,
                "pred": early_v_mid_late_pred,
                "actual": early_v_mid_late_actual,
            }
        )
        
        df_pred_out_mean = df_pred_out.groupby('accession_id').mean().reset_index()
        df_pred_out.to_csv(results_dir / f"holdout_predictions_with_acc_ids_top_n_sig_{top_n_sig}.csv", index=False)
        df_pred_out_mean.to_csv(results_dir / f"holdout_predictions_agg_with_acc_ids_top_n_sig_{top_n_sig}.csv", index=False)
        fpr, tpr, thresholds = roc_curve(early_v_mid_late_actual, early_v_mid_late_pred)
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, color="blue", label="ROC curve")
        plt.plot([0, 1], [0, 1], color="red", linestyle="--", label="Random guess")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate (Recall)")
        plt.title("ROC Curve")
        plt.legend()
        plt.savefig(results_dir / f"holdout_roc_curve_top_n_sig_{top_n_sig}.png")
        print('THRESHOLDS')
        print(thresholds)
        ras = roc_auc_score(early_v_mid_late_actual, early_v_mid_late_pred)

        print(f"ROC AUC score ({embeddings_in}): {ras}")

        fpr, tpr, thresholds = roc_curve(df_pred_out_mean.actual, df_pred_out_mean.pred)
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, color="blue", label="ROC curve")
        plt.plot([0, 1], [0, 1], color="red", linestyle="--", label="Random guess")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate (Recall)")
        plt.title("ROC Curve")
        plt.legend()
        plt.savefig(results_dir / f"holdout_agg_roc_curve_top_n_sig_{top_n_sig}.png")
        ras = roc_auc_score(df_pred_out_mean.actual, df_pred_out_mean.pred)

        print(f"ROC AUC score ({embeddings_in}): {ras}")


def main(embedding_folder):
    results_dir = Path(f"modeling/{str(datetime.now()).replace(' ','_')}")
    results_dir.mkdir(exist_ok=True)
    embeddings_train = {}
    embeddings_holdout = {}
    prepare_train = {}
    prepare_holdout = {}
    top_n_sig = None
    
    for file_name in embedding_folder.iterdir():
        if file_name.is_file():
            if file_name.suffix == ".json":
                print(file_name)
                with open(file_name, "r") as fp:
                    metadata = json.load(fp)
                kfold = int(file_name.stem.split('kfold_')[-1].split('_')[0])
                
                if embeddings_train == {}:
                    for top_n_sig in metadata['top_n_sig']:
                        embeddings_train[top_n_sig] = {}
                        embeddings_holdout[top_n_sig] = {}
                        prepare_train[top_n_sig] = {}
                        prepare_holdout[top_n_sig] = {}
                print(f'Kfold = {kfold}')

                for top_n_sig, train_embedding_filename, holdout_embedding_filename, train_prepare_filename, holdout_prepare_filename in zip(
                    metadata['top_n_sig'], 
                    metadata['train_embedding_filenames'], 
                    metadata['holdout_embedding_filenames'], 
                    metadata['train_annotated_filenames'], 
                    metadata['holdout_annotated_filenames']):
                    kfold = int(train_embedding_filename.split('kfold_')[-1].split('_')[0])
                   
                    embeddings_train[top_n_sig][kfold] = train_embedding_filename
                    embeddings_holdout[top_n_sig][kfold] = holdout_embedding_filename
                    prepare_train[top_n_sig][kfold] = train_prepare_filename
                    prepare_holdout[top_n_sig][kfold] = holdout_prepare_filename
                '''
                for i, top_n_sig in enumerate(metadata["top_n_sig"]):
                    if top_n_sig in embeddings.keys():
                        
                    else:
                        if 'train' in embeddings[top_n_sig].keys():
                                        
                embeddings_filenames = []
                for embedding_train_fname in metadata["embedding_train_fnames"]:
                    print(f"Running modeling on {embedding_train_fname}")
                    embeddings_filename = c2s(embedding_train_fname, results_dir)
                    embeddings_filenames.append(str(embeddings_filename))

                metadata["modeling_filenames"] = str(embeddings_filenames)

                with open(results_dir / f"modeling_{file_name.stem}.json", "w") as fp:
                    json.dump(metadata, fp)
                '''
                deg_results_path = Path(metadata['deg_results_file'])
    with open(deg_results_path.parent / 'fold_indices.pkl', 'rb') as file:
        # Load and reconstruct the object
        kfolds_in = pickle.load(file)
    for top_n_sig in embeddings_train.keys():
        print(f'modeling for top n {top_n_sig}')
        model(top_n_sig, embeddings_train[top_n_sig], prepare_train[top_n_sig], results_dir, kfolds_in=kfolds_in, holdout_embeddings_in=embeddings_holdout[top_n_sig], prepare_holdout = prepare_holdout[top_n_sig])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_modeling.py <embedding_path_in>")
        sys.exit(1)
    main(Path(sys.argv[1]))