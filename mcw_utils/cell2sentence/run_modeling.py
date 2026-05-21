import sys
from sklearn.manifold import TSNE
import plotly.express as px
import pandas as pd

import numpy as np
from FCNet_prekfold import train_fc_network_prekfold
from FCNet import train_fc_network
from FCNet_single import train_fc_network_single
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt


def embedding_to_X_y(embeddings_in, max_layer):

    embeddings_df = pd.read_csv(f"embeddings/{embeddings_in}.csv")
    embeddings_df.set_index(["accession_id", "layer_id"], inplace=True)
    label_data = pd.read_csv("/mnt/c/Users/msochor/Downloads/big_data_with_labels.csv")

    early_mid_late_data = label_data[
        label_data.recurrence_time_sur.isin(["early", "mid", "late"])
    ]
    early_mid_late_panc_data = early_mid_late_data[early_mid_late_data.is_panc == True]

    tsne = TSNE(n_components=2, random_state=0, perplexity=5)

    proj = tsne.fit_transform(embeddings_df)
    fig = px.scatter(
        proj,
        x=0,
        y=1,
        color=embeddings_df.reset_index()["layer_id"],
        labels={"color": "layer_id"},
    )

    fig.write_html(f"tsne_plots/tsne_plot_{embeddings_in}.html")

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


def main(embeddings_in, kfolds_in=None, holdout_embeddings_in=None):
    max_layer = 9
    if kfolds_in is not None:
        Xs = []
        ys = []
        acc_ids = []
        for embedding_in in embeddings_in.split(","):
            X, y, acc_id = embedding_to_X_y(embedding_in, max_layer)
            Xs.append(X)
            ys.append(y)
            acc_ids.append(acc_id)
        models, histories, all_preds, all_preds_proba, all_y_vals, all_acc_ids = (
            train_fc_network_prekfold(
                Xs,
                ys,
                acc_ids,
                kfolds_in,
                hidden_sizes=(512, 128),
                n_classes=3,
                epochs=15,
                batch_size=2,
                lr=1e-3,
            )
        )
    else:
        X, y, acc_ids = embedding_to_X_y(embeddings_in, max_layer)
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
    df_pred_out.to_csv(f"predictions_with_acc_ids_{embeddings_in}.csv", index=False)
    fpr, tpr, thresholds = roc_curve(early_v_mid_late_actual, early_v_mid_late_pred)
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color="blue", label="ROC curve")
    plt.plot([0, 1], [0, 1], color="red", linestyle="--", label="Random guess")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate (Recall)")
    plt.title("ROC Curve")
    plt.legend()
    plt.savefig(f"roc_curve_{embeddings_in}.png")
    ras = roc_auc_score(early_v_mid_late_actual, early_v_mid_late_pred)

    print(f"ROC AUC score ({embeddings_in}): {ras}")
    if holdout_embeddings_in is not None:

        holdout_embeddings_df = pd.read_csv(f"{holdout_embeddings_in}.csv")
        holdout_embeddings_df.set_index(["accession_id", "layer_id"], inplace=True)
        label_data_holdout = pd.read_csv(
            "/mnt/c/Users/msochor/Downloads/big_data_with_labels_holdout.csv"
        )

        early_mid_late_data_holdout = label_data_holdout[
            label_data_holdout.recurrence_time_sur.isin(["early", "mid", "late"])
        ]
        early_mid_late_panc_data_holdout = early_mid_late_data_holdout[
            early_mid_late_data_holdout.is_panc == True
        ]
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

        val = history["val_acc"][-1]
        print(f"Max layers: {max_layer}")
        print(f"Average final validation acc: {val}")

        print(f"Confusion matrix for holdout:")
        cm = confusion_matrix(y_vals, preds)
        print(cm)

        early_v_mid_late_actual = []
        for y_val in y_vals:
            if y_val == 0:
                early_v_mid_late_actual.append(0)
            else:
                early_v_mid_late_actual.append(1)
        early_v_mid_late_pred = []
        for y_pred in preds_proba:
            early_v_mid_late_pred.append(1 - y_pred[0])
        fpr, tpr, thresholds = roc_curve(early_v_mid_late_actual, early_v_mid_late_pred)
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, color="blue", label="ROC curve")
        plt.plot([0, 1], [0, 1], color="red", linestyle="--", label="Random guess")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate (Recall)")
        plt.title("ROC Curve")
        plt.legend()
        plt.savefig(f"roc_curve_{holdout_embeddings_in}_holdout.png")
        ras = roc_auc_score(early_v_mid_late_actual, early_v_mid_late_pred)

        print(f"ROC AUC score ({holdout_embeddings_in}): {ras}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(
            "Usage: python run_modeling.py <embeddings_in> (optional)<kfolds_in> (optional) holdout_embeddings_in"
        )
        sys.exit(1)
    if len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[1], kfolds_in=sys.argv[2])
    else:
        if sys.argv[2].lower() == "none":
            kfolds_in = None
        else:
            kfolds_in = sys.argv[2].lower()
        main(sys.argv[1], kfolds_in=kfolds_in, holdout_embeddings_in=sys.argv[3])
