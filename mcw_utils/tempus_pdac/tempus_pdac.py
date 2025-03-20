import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import pandas as pd
from ffile import Ffile

from ._util import trim_emrn, trim_hashed_patient_name, trim_mrn, make_tall_variants


@dataclass
class TempusPdac:
    pdac_sotb_fname: str
    tempus_dirs: list[str]
    tempus_manifest_fname: pd.DataFrame
    rcc_fname: str
    grey: str = "rgba(30, 30, 30, 0.2)"
    blue: str = "rgba(0, 0, 180, 0.4)"
    red: str = "rgba(180, 0, 0, 0.4)"
    yellow: str = "rgba(180, 180, 0, 0.4)"
    green: str = "rgba(0, 180, 0, 0.4)"

    from ._qc import (
        _find_missing_files,
        _quantify_patient_level_merge,
        _quantify_rcc_merge,
        _quantify_report_level_merge,
    )
    from ._read import _read_pdac_sotb, _read_rcc, _read_tempus, _read_tempus_manifest

    def __post_init__(self):
        now = datetime.now()

        # Extract year, month, and day
        year = now.year
        month = now.month
        day = now.day
        self.outdir = Path(f"./output_{month}_{day}_{year}/")
        self.outdir.mkdir(exist_ok=True)
        input_files = {
            "pdac_sotb_fname": self.pdac_sotb_fname,
            "tempus_dirs": self.tempus_dirs,
            "tempus_manifest_fname": self.tempus_manifest_fname,
            "rcc_fname": self.rcc_fname,
        }
        self.input_files_str = json.dumps(input_files)
        with open(self.outdir / "input_files.json", "w") as fp:
            fp.write(self.input_files_str)

    def read_files(self):
        self._read_pdac_sotb()
        self._read_tempus()
        self._read_tempus_manifest()
        self._read_rcc()

    def merge(self):
        tempus_and_pdac = self.tempus_meta.merge(
            self.pdac_sotb_min, how="left", left_on="emr_id", right_on="emrn"
        )
        tempus_and_pdac_2 = tempus_and_pdac.merge(
            self.pdac_sotb_min, how="left", left_on="integer_emr_id", right_on="mrn"
        )

        tempus_and_pdac_2["mrn"] = tempus_and_pdac_2.apply(trim_mrn, axis=1)
        tempus_and_pdac_2["emrn"] = tempus_and_pdac_2.apply(trim_emrn, axis=1)
        tempus_and_pdac_2["hashed_patient_name"] = tempus_and_pdac_2.apply(
            trim_hashed_patient_name, axis=1
        )

        tempus_and_pdac_2.drop(
            columns=[
                "mrn_x",
                "emrn_x",
                "missing_mrn_x",
                "missing_emrn_x",
                "mrn_y",
                "emrn_y",
                "missing_mrn_y",
                "missing_emrn_y",
                "hashed_patient_name_x",
                "hashed_patient_name_y",
            ],
            inplace=True,
        )
        self.tempus_and_pdac = tempus_and_pdac_2

        self.tempus_and_pdac_full = self.pdac_sotb_min.merge(
            tempus_and_pdac_2, on="mrn", how="left"
        )
        self.tempus_pdac_manifest = tempus_and_pdac_2.merge(
            self.tempus_manifest,
            how="left",
            left_on="accession_id",
            right_on="Accession Number",
        )
        self.tempus_pdac_manifest_rcc = self.tempus_pdac_manifest.merge(
            self.rcc,
            how="left",
            left_on=["Accession Number", "report_type"],
            right_on=["accession_number", "report_type"],
        )

        self.tempus_pdac_manifest_rcc.fillna(
            value={"has_rcc_bam_file": False}, inplace=True
        )
        self.tempus_pdac_manifest_rcc.to_csv(
            self.outdir / "tempus_json_pdac_rcc_merged.csv", index=False
        )
        self.tempus_pdac_manifest_rcc_dna_tall = make_tall_variants(
            self.tempus_pdac_manifest_rcc
        )
        self.tempus_pdac_manifest_rcc_dna_tall.to_csv(
            self.outdir / "tempus_json_pdac_rcc_merged_dna_tall.csv", index=False
        )
        tempus_pdac_with_meta_has_rcc = self.tempus_pdac_manifest_rcc[
            self.tempus_pdac_manifest_rcc.has_rcc_bam_file == True
        ]
        rna_tempus_pdac_with_meta_has_rcc = tempus_pdac_with_meta_has_rcc[
            tempus_pdac_with_meta_has_rcc.report_type == "RNA"
        ]
        rna_tempus_pdac_with_meta_has_rcc = rna_tempus_pdac_with_meta_has_rcc[
            rna_tempus_pdac_with_meta_has_rcc.mrn.notnull()
        ]

        rna_tempus_pdac_with_meta_has_rcc.to_csv(
            self.outdir / "RNA_rcc_bam_with_pdac_sotb_mrn.csv", index=False
        )
        # run merge quality controls
        self._find_missing_files()
        self._quantify_patient_level_merge()
        self._quantify_report_level_merge()
        self._quantify_rcc_merge()

        local_dir = Path(__file__).parent
        report_template = Ffile(local_dir / "report.html.template")
        with open(self.outdir / "results.html", "w") as fp:
            fp.write(
                report_template.f(
                    pdac_sotb_fname=self.pdac_sotb_fname,
                    tempus_dirs="<br>".join(self.tempus_dirs),
                    tempus_manifest_fname=self.tempus_manifest_fname,
                    rcc_fname=self.rcc_fname,
                    fig1=self.fig1,
                    fig2=self.fig2,
                    fig3=self.fig3,
                    fig4=self.fig4,
                )
            )
