import sys
from sklearn.manifold import TSNE
import plotly.express as px
import pandas as pd

import numpy as np
from FCNet import train_fc_network
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt


def main(embeddings_in):
    embeddings_df = pd.read_csv(f"{embeddings_in}.csv")
    embeddings_df.set_index(["accession_id", "layer_id"], inplace=True)
    label_data = pd.read_csv("/mnt/c/Users/msochor/Downloads/big_data_with_labels.csv")
    model_data = pd.read_csv(
        "/mnt/c/Users/msochor/Downloads/big_data_for_c2s_modeling.csv"
    )

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

    fig.write_html(f"tsne_plot_{embeddings_in}.html")

    max_layer = 9
    X = []
    y = []
    for i in range(len(early_mid_late_panc_data.accession_id.values)):
        row = early_mid_late_panc_data.iloc[i]
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
    models, histories, all_preds, all_preds_proba, all_y_vals = train_fc_network(
        X,
        y,
        hidden_sizes=(512, 128),
        n_classes=3,
        epochs=15,
        batch_size=2,
        lr=1e-3,
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


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_modeling.py <embeddings_in>")
        sys.exit(1)
    main(sys.argv[1])
