import math
from hashlib import sha256

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plottable import Table
from pathlib import Path


def trim_mrn(row):
    if np.isnan(row.mrn_x):
        if np.isnan(row.mrn_y):
            return row.mrn_y
        else:
            return int(row.mrn_y)
    else:
        if np.isnan(row.mrn_x):
            return row.mrn_x
        else:
            return int(row.mrn_x)


def has_json_file(row):
    if isinstance(row.acc_num, float) and math.isnan(row.acc_num):
        return False
    else:
        return True


def hash_patient_name(row):
    return sha256(bytes(row.first_name + row.last_name, "utf-8")).hexdigest()


def emr_id_to_integer(row):
    if row.emr_id[0] == "E":
        return int(row.emr_id[1:])
    else:
        return int(row.emr_id)


def trim_emrn(row):
    if np.isnan(row.mrn_x):
        return row.emrn_y
    else:
        return row.emrn_x


def trim_hashed_patient_name(row):
    if isinstance(row.hashed_patient_name_x, float) and math.isnan(
        row.hashed_patient_name_x
    ):
        return row.hashed_patient_name_y
    else:
        return row.hashed_patient_name_x


def plot_sankey(title: str, nodes: dict, connections: dict, outdir: Path = None):
    labels = []
    label_colors = []
    for node in nodes:
        labels.append(node[0])
        label_colors.append(node[1])

    sources = []
    targets = []
    values = []
    connection_colors = []
    for connection in connections:
        source, target = connection[0].split("->")
        sources.append(source)
        targets.append(target)
        values.append(connection[1])
        connection_colors.append(connection[2])

    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.5),
                    label=labels,
                    color=label_colors,
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=connection_colors,
                ),
            )
        ]
    )

    fig.update_layout(title_text=title, font_size=10)
    outfile = title.lower().replace(" ", "_") + "_sankey.png"
    if outdir:
        outfile = outdir / outfile

    fig.show()
    fig.write_image(outfile)
    return fig.to_html(full_html=False)


def plot_table(
    title,
    count_data: list[list[int]],
    column_names: list[str],
    row_names: list[str],
    outdir: Path = None,
):

    total_count_data = []
    for count_data_row in count_data:
        total_count_data.append(count_data_row + [sum(count_data_row)])
    # data = np.array(total_count_data)

    df = pd.DataFrame(
        total_count_data, columns=column_names + ["Total"], index=row_names
    )
    df.index.name = ""
    # Create the figure and axes
    fig, ax = plt.subplots()

    table = Table(
        df,
        textprops={"ha": "center"},
        row_dividers=False,  # remove lines between rows
        col_label_divider=True,
        footer_divider=True,
        odd_row_color="lightgrey",
    )
    table.columns["Total"].set_facecolor("lightblue")
    outfile = title.lower().replace(" ", "_") + "_table.png"
    if outdir:
        outfile = outdir / outfile
    plt.savefig(outfile)
    plt.show()


def count_dna_reports(df):
    return len(df[df.report_type == "DNA"])


def count_rna_reports(df):
    return len(df[df.report_type == "RNA"])


def make_tall_variants(df):
    c = []
    for col in df.columns:
        if col != "variants":
            c.append(col)
    tall = []
    for i in range(len(df)):
        if isinstance(df.iloc[i].variants, str) and not df.iloc[i].variants == "":
            variants = df.iloc[i].variants.split("|")
            for variant in variants:
                gene, mut, source, comment, allelic_fraction = variant.split(":")
                tall.append(
                    pd.concat(
                        [
                            df.iloc[i][c],
                            pd.Series(
                                {
                                    "gene": gene,
                                    "variant": mut,
                                    "source": source,
                                    "comment": comment,
                                    "allelic_fraction": allelic_fraction,
                                }
                            ),
                        ]
                    )
                )
    return pd.concat(tall, axis=1).T
