import numpy as np
from mcw_utils.tempus_pdac._util import (
    count_dna_reports,
    count_rna_reports,
    has_json_file,
    plot_sankey,
    plot_table,
)


def _find_missing_files(self):
    self.missing_mrn = self.tempus_and_pdac[self.tempus_and_pdac.mrn.isnull()]
    missing_mrn_unique = (
        self.missing_mrn.groupby("hashed_patient_id").first().reset_index()
    )
    missing_mrn_unique = (
        self.missing_mrn.groupby("hashed_patient_id").first().reset_index()
    )

    # Find missing RNA and DNA files and save to csv

    has_rna_list = self.tempus_pdac_manifest[
        self.tempus_pdac_manifest["Has RNA Folder"] == "Yes"
    ]
    has_rna_list_remove_dna = has_rna_list[has_rna_list.report_type == "RNA"]
    accession_id_missing_rna_json = set(has_rna_list.accession_id) - set(
        has_rna_list_remove_dna.accession_id
    )
    self.tempus_pdac_missing_rna_json = self.tempus_pdac_manifest[
        self.tempus_pdac_manifest.accession_id.isin(list(accession_id_missing_rna_json))
    ]
    self.missing_rna_json_count = len(accession_id_missing_rna_json)

    has_dna_list = self.tempus_pdac_manifest[
        self.tempus_pdac_manifest["Has DNA Folder"] == "Yes"
    ]
    has_dna_list_remove_rna = has_dna_list[has_dna_list.report_type == "DNA"]
    accession_id_missing_dna_json = set(has_dna_list.accession_id) - set(
        has_dna_list_remove_rna.accession_id
    )
    self.missing_dna_json_count = len(accession_id_missing_dna_json)

    self.tempus_pdac_missing_dna_json = self.tempus_pdac_manifest[
        self.tempus_pdac_manifest.accession_id.isin(list(accession_id_missing_dna_json))
    ]

    missing_mrn_unique.to_csv(
        self.outdir / "patients_with_json_missing_from_pdac_data_pull.csv"
    )
    self.tempus_pdac_missing_dna_json.to_csv(
        self.outdir / "tempus_pdac_merged_missing_dna_json.csv"
    )
    self.tempus_pdac_missing_rna_json.to_csv(
        self.outdir / "tempus_pdac_merged_missing_rna_json.csv"
    )


def _quantify_report_level_merge(self):
    dna_report_count = len(
        self.tempus_pdac_manifest[self.tempus_pdac_manifest.report_type == "DNA"]
    )
    rna_report_count = len(
        self.tempus_pdac_manifest[self.tempus_pdac_manifest.report_type == "RNA"]
    )

    dna_report_with_rna_count = len(
        self.tempus_pdac_manifest[
            np.logical_and(
                self.tempus_pdac_manifest.report_type == "DNA",
                self.tempus_pdac_manifest["Has RNA Folder"] == "Yes",
            )
        ]
    )

    dna_report_without_rna_count = len(
        self.tempus_pdac_manifest[
            np.logical_and(
                self.tempus_pdac_manifest.report_type == "DNA",
                self.tempus_pdac_manifest["Has RNA Folder"] == "No",
            )
        ]
    )

    rna_report_with_dna_count = len(
        self.tempus_pdac_manifest[
            np.logical_and(
                self.tempus_pdac_manifest.report_type == "RNA",
                self.tempus_pdac_manifest["Has DNA Folder"] == "Yes",
            )
        ]
    )

    rna_report_without_dna_count = len(
        self.tempus_pdac_manifest[
            np.logical_and(
                self.tempus_pdac_manifest.report_type == "RNA",
                self.tempus_pdac_manifest["Has DNA Folder"] == "No",
            )
        ]
    )

    missing_dna_report_with_rna_count = len(
        self.tempus_pdac_missing_dna_json[
            self.tempus_pdac_missing_dna_json["Has RNA Folder"] == "Yes"
        ]
    )
    missing_dna_report_without_rna_count = len(
        self.tempus_pdac_missing_dna_json[
            self.tempus_pdac_missing_dna_json["Has RNA Folder"] == "No"
        ]
    )

    missing_rna_report_with_dna_count = len(
        self.tempus_pdac_missing_rna_json[
            self.tempus_pdac_missing_rna_json["Has DNA Folder"] == "Yes"
        ]
    )
    missing_rna_report_without_dna_count = len(
        self.tempus_pdac_missing_rna_json[
            self.tempus_pdac_missing_rna_json["Has DNA Folder"] == "No"
        ]
    )
    nodes = [
        ["Expected Tempus JSON from manifest", "grey"],
        ["Missing Tempus JSON - DNA", "red"],
        ["Missing Tempus JSON - RNA", "red"],
        ["Actual Tempus JSON - DNA", "blue"],
        ["Actual Tempus JSON - RNA", "yellow"],
        ["Expected Tempus JSON - Patient DNA and RNA", "green"],
        ["Expected Tempus JSON - Patient DNA only", "blue"],
        ["Expected Tempus JSON - Patient RNA only", "yellow"],
    ]
    connections = [
        ["0->1", self.missing_dna_json_count, self.red],
        ["0->2", self.missing_rna_json_count, self.red],
        ["0->3", dna_report_count, self.blue],
        ["0->4", rna_report_count, self.yellow],
        ["3->5", dna_report_with_rna_count, self.blue],
        ["4->5", rna_report_with_dna_count, self.yellow],
        ["3->6", dna_report_without_rna_count, self.blue],
        ["4->7", rna_report_without_dna_count, self.yellow],
        ["1->5", missing_dna_report_with_rna_count, self.red],
        ["1->6", missing_dna_report_without_rna_count, self.red],
        ["2->5", missing_rna_report_with_dna_count, self.red],
        ["2->7", missing_rna_report_without_dna_count, self.red],
    ]
    title = "Reports with Tempus files"
    self.fig3 = plot_sankey(title, nodes, connections, outdir=self.outdir)
    data = [
        [
            dna_report_count,
            rna_report_count,
            self.missing_dna_json_count,
            self.missing_rna_json_count,
        ],
        [
            dna_report_with_rna_count,
            rna_report_with_dna_count,
            missing_dna_report_with_rna_count,
            missing_rna_report_with_dna_count,
        ],
        [
            dna_report_without_rna_count,
            0,
            missing_dna_report_without_rna_count,
            0,
        ],
        [
            0,
            rna_report_without_dna_count,
            0,
            missing_rna_report_without_dna_count,
        ],
    ]

    columns = ["DNA files", "RNA files", "Missing\nDNA", "Missing\nRNA"]
    rows = [
        "Expected\nJSONs from\nmanifest",
        "Patients\nwith DNA\nand RNA",
        "Patients\nwith only\nDNA",
        "Patients\nwith only\nRNA",
    ]
    plot_table(title, data, columns, rows, outdir=self.outdir)


def _quantify_patient_level_merge(self):

    patients_in_pdac_missing_tempus = len(
        set(
            self.tempus_and_pdac_full[
                self.tempus_and_pdac_full.emr_id.isnull()
            ].hashed_patient_name_x
        )
    )
    patients_in_tempus_and_pdac = len(
        set(
            self.tempus_and_pdac[
                self.tempus_and_pdac.hashed_patient_name.notnull()
            ].hashed_patient_name
        )
    )
    patients_in_tempus_missing_pdac = len(set(self.missing_mrn.hashed_patient_id))

    dna_report_count = (
        self.tempus_pdac_manifest.groupby("hashed_patient_id")
        .apply(count_dna_reports)
        .reset_index()
        .rename(columns={0: "dna_report_count"})
    )
    rna_report_count = (
        self.tempus_pdac_manifest.groupby("hashed_patient_id")
        .apply(count_rna_reports)
        .reset_index()
        .rename(columns={0: "rna_report_count"})
    )
    report_counts = dna_report_count.merge(rna_report_count, on="hashed_patient_id")
    patients_with_dna_reports_only = len(
        report_counts[
            np.logical_and(
                report_counts.rna_report_count == 0,
                report_counts.dna_report_count > 0,
            )
        ]
    )
    patients_with_rna_reports_only = len(
        report_counts[
            np.logical_and(
                report_counts.dna_report_count == 0,
                report_counts.rna_report_count > 0,
            )
        ]
    )
    patients_with_dna_and_rna_reports = len(
        report_counts[
            np.logical_and(
                report_counts.dna_report_count > 0,
                report_counts.rna_report_count > 0,
            )
        ]
    )
    nodes = [
        ["PDAC SOTB data pull", "grey"],
        ["Missing from SOTB pull", "red"],
        ["PDAC SOTB patients - no Tempus", "grey"],
        ["PDAC SOTB patients - Tempus JSON", "green"],
    ]

    connections = [
        ["0->2", patients_in_pdac_missing_tempus, self.grey],
        ["0->3", patients_in_tempus_and_pdac, self.blue],
        ["1->3", patients_in_tempus_missing_pdac, self.red],
    ]
    title = "Patient count PDAC tempus"
    self.fig1 = plot_sankey(title, nodes, connections, outdir=self.outdir)

    count_data = [
        [
            patients_in_tempus_and_pdac,
            patients_in_pdac_missing_tempus,
            patients_in_tempus_missing_pdac,
        ],
    ]
    columns = ["PDAC and\n Tempus", "PDAC only", "Tempus only"]
    rows = ["Patient count"]
    plot_table(title, count_data, columns, rows, outdir=self.outdir)

    nodes = [
        ["PDAC SOTB patients - Tempus JSON", "green"],
        ["Patients with DNA Json only", "blue"],
        ["Patients with RNA Json only", "yellow"],
        ["Patients with DNA and RNA Json", "green"],
    ]

    connections = [
        ["0->1", patients_with_dna_reports_only, self.blue],
        ["0->2", patients_with_rna_reports_only, self.yellow],
        ["0->3", patients_with_dna_and_rna_reports, self.green],
    ]
    title = "Tempus Patient count DNA vs RNA"
    self.fig2 = plot_sankey(title, nodes, connections, outdir=self.outdir)
    count_data = [
        [
            patients_with_dna_reports_only,
            patients_with_rna_reports_only,
            patients_with_dna_and_rna_reports,
        ]
    ]
    columns = ["DNA only", "RNA only", "DNA and RNA"]
    rows = ["Tempus Patient\n count"]
    plot_table(title, count_data, columns, rows, outdir=self.outdir)


def _quantify_rcc_merge(self):

    tempus_pdac_manifest_missing_rcc = self.tempus_pdac_manifest_rcc[
        self.tempus_pdac_manifest_rcc.has_rcc_bam_file == False
    ]

    count_tempus_dna_no_rcc = len(
        tempus_pdac_manifest_missing_rcc[
            tempus_pdac_manifest_missing_rcc.report_type == "DNA"
        ]
    )
    count_tempus_rna_no_rcc = len(
        tempus_pdac_manifest_missing_rcc[
            tempus_pdac_manifest_missing_rcc.report_type == "RNA"
        ]
    )

    tempus_pdac_manifest_missing_rcc.to_csv(
        self.outdir / "rcc_bam_files_missing_tempus_json.csv"
    )

    rcc_tempus_pdac_with_manifest = self.rcc.merge(
        self.tempus_pdac_manifest,
        how="left",
        right_on=["Accession Number", "report_type"],
        left_on=["accession_number", "report_type"],
    )

    rcc_tempus_pdac_with_manifest["has_json_file"] = (
        rcc_tempus_pdac_with_manifest.apply(has_json_file, axis=1)
    )
    tempus_pdac_manifest_rcc_missing_json = rcc_tempus_pdac_with_manifest[
        rcc_tempus_pdac_with_manifest.has_json_file == False
    ]
    tempus_pdac_manifest_rcc_has_json = rcc_tempus_pdac_with_manifest[
        rcc_tempus_pdac_with_manifest.has_json_file == True
    ]

    count_rcc_missing_json_rna = len(
        tempus_pdac_manifest_rcc_missing_json[
            tempus_pdac_manifest_rcc_missing_json.report_type == "RNA"
        ]
    )
    count_rcc_missing_json_dna = len(
        tempus_pdac_manifest_rcc_missing_json[
            tempus_pdac_manifest_rcc_missing_json.report_type == "DNA"
        ]
    )
    tempus_pdac_manifest_rcc_missing_json.to_csv(
        self.outdir / "tempus_json_files_missing_rcc_bam.csv"
    )

    count_rcc_has_json_rna = len(
        tempus_pdac_manifest_rcc_has_json[
            tempus_pdac_manifest_rcc_has_json.report_type == "RNA"
        ]
    )
    count_rcc_has_json_dna = len(
        tempus_pdac_manifest_rcc_has_json[
            tempus_pdac_manifest_rcc_has_json.report_type == "DNA"
        ]
    )
    nodes = [
        ["BAM files on RCC", "green"],
        ["RCC only - missing Tempus JSON DNA", "red"],
        ["RCC only - missing Tempus JSON RNA", "red"],
        ["RCC + Tempus JSON - DNA", "blue"],
        ["RCC + Tempus JSON - RNA", "yellow"],
        ["Total RCC + Tempus JSON", "green"],
        ["Total RCC - missing Tempus JSON", "red"],
        ["Tempus JSON - missing RCC DNA", "red"],
        ["Tempus JSON - missing RCC RNA", "red"],
        ["Total Tempus JSON", "green"],
    ]
    connections = [
        ["0->1", count_rcc_missing_json_dna, self.red],
        ["0->2", count_rcc_missing_json_rna, self.red],
        ["0->3", count_rcc_has_json_dna, self.blue],
        ["0->4", count_rcc_has_json_rna, self.yellow],
        ["3->5", count_rcc_has_json_dna, self.blue],
        ["4->5", count_rcc_has_json_rna, self.yellow],
        ["1->6", count_rcc_missing_json_dna, self.red],
        ["2->6", count_rcc_missing_json_rna, self.red],
        ["5->9", count_rcc_has_json_dna + count_rcc_has_json_rna, self.green],
        ["7->9", count_tempus_dna_no_rcc, self.red],
        ["8->9", count_tempus_rna_no_rcc, self.red],
    ]
    title = "Tempus JSON joining to RCC BAM files"
    self.fig4 = plot_sankey(title, nodes, connections, outdir=self.outdir)
    data = [
        [
            count_rcc_has_json_dna,
            count_rcc_has_json_rna,
            count_tempus_dna_no_rcc,
            count_tempus_rna_no_rcc,
        ],
        [
            count_rcc_has_json_dna,
            count_rcc_has_json_rna,
            count_rcc_missing_json_dna,
            count_rcc_missing_json_rna,
        ],
    ]
    columns = ["DNA files", "RNA files", "Missing\nDNA", "Missing\nRNA"]
    rows = ["Tempus JSONS\ncount", "RCC BAM\nfile count"]
    plot_table(title, data, columns, rows, outdir=self.outdir)
