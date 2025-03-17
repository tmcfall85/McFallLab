import json
from hashlib import sha256
from pathlib import Path

import pandas as pd

from mcw_utils.tempus_pdac._util import emr_id_to_integer, hash_patient_name


def _read_tempus_files_in_directory(directory_path):
    """Reads and returns the content of all files in a directory.

    Args:
        directory_path: Path to the directory (string or Path object).

    Returns:
        A dictionary where keys are file names and values are file contents.
    """
    directory = Path(directory_path)
    if not directory.is_dir():
        raise NotADirectoryError(f"{directory_path} is not a directory")
    emr_id = []
    hashed_patient_id = []
    schema_version = []
    acc_num = []
    accession_id = []
    report_id = []
    kras_variants = []
    report_type = []
    bio_inf_pipeline_version = []
    tempus_id = []
    test_code = []
    specimen_count = []
    specimen_sample_category = []
    specimen_sample_site = []
    specimen_block_id = []
    specimen_tumor_percentage = []

    for file_path in directory.iterdir():
        if file_path.is_dir():
            for inner_file_path in file_path.iterdir():
                if inner_file_path.is_file() and inner_file_path.suffix == ".json":
                    has_json = True
                    with open(inner_file_path, "r") as file:
                        data = json.loads(file.read())
                        if data["metadata"]["schemaVersion"] in ["1.3.1", "1.3.2"]:
                            emr_id.append(data["patient"]["emr_id"])
                        else:
                            emr_id.append(data["patient"]["emrId"])

                        hashed_patient_id.append(
                            sha256(
                                bytes(
                                    data["patient"]["firstName"]
                                    + data["patient"]["lastName"],
                                    "utf-8",
                                )
                            ).hexdigest()
                        )
                        # hashed_patient_id.append(data['patient']['firstName'] + data['patient']['lastName'])
                        schema_version.append(data["metadata"]["schemaVersion"])
                        acc_num.append(file_path.parts[-1])
                        accession_id.append(data["order"]["accessionId"])
                        report_type.append(data["report"]["workflow"]["reportType"])
                        has_kras_variant = False
                        if report_type == "DNA":
                            muts = data["results"][
                                "somaticPotentiallyActionableMutations"
                            ]
                            for mut in muts:
                                if mut["gene"] == "KRAS":
                                    mut_kras_variants = []
                                    for variant in mut["gene"]["variants"]:
                                        mut_kras_variants.append(
                                            variant["mutationEffect"]
                                        )
                                    has_kras_variant = True
                                    kras_variants.append("_".join(mut_kras_variants))
                        if has_kras_variant == False:
                            kras_variants.append("")

                        report_id.append(data["report"]["reportId"])
                        bio_inf_pipeline_version.append(
                            data["report"]["bioInfPipeline"]
                        )
                        tempus_id.append(data["patient"]["tempusId"])
                        test_code.append(data["order"]["test"]["code"])
                        specimen_count.append(len(data["specimens"]))
                        specimen_sample_category.append(
                            data["specimens"][0]["sampleCategory"]
                        )
                        specimen_sample_site.append(data["specimens"][0]["sampleSite"])
                        specimen_block_id.append(
                            data["specimens"][0]["institutionData"]["blockId"]
                        )
                        specimen_tumor_percentage.append(
                            data["specimens"][0]["institutionData"]["tumorPercentage"]
                        )

    df = pd.DataFrame(
        {
            "emr_id": emr_id,
            "hashed_patient_id": hashed_patient_id,
            "schema_version": schema_version,
            "acc_num": acc_num,
            "accession_id": accession_id,
            "report_type": report_type,
            "report_id": report_id,
            "bio_inf_pipeline_version": bio_inf_pipeline_version,
            "tempus_id": tempus_id,
            "test_code": test_code,
            "kras_variants": kras_variants,
            "specimen_count": specimen_count,
            "specimen_sample_category": specimen_sample_category,
            "specimen_sample_site": specimen_sample_site,
            "specimen_block_id": specimen_block_id,
            "specimen_tumor_percentage": specimen_tumor_percentage,
        }
    )

    return df


def _read_pdac_sotb(self):
    pdac_sotb = pd.read_stata(self.pdac_sotb_fname)
    pdac_sotb["hashed_patient_name"] = pdac_sotb.apply(hash_patient_name, axis=1)
    pdac_sotb_min = pdac_sotb[["mrn", "emrn", "hashed_patient_name"]].copy()
    pdac_sotb_min["missing_mrn"] = pdac_sotb_min.mrn == 1
    pdac_sotb_min["missing_emrn"] = pdac_sotb_min.emrn == ""

    self.pdac_sotb_min = pdac_sotb_min


def _read_tempus(self):
    meta_data = []
    for tempus_dir in self.tempus_dirs:
        meta_data.append(_read_tempus_files_in_directory(tempus_dir))

    tempus_meta = pd.concat(meta_data)
    tempus_meta["integer_emr_id"] = tempus_meta.apply(emr_id_to_integer, axis=1)
    self.tempus_meta = tempus_meta


def _read_tempus_manifest(self):
    self.tempus_manifest = pd.read_csv(self.tempus_manifest_fname)


def _read_rcc(self):
    with open(self.rcc_fname, "r") as fp:
        rcc_files = fp.readlines()

    accession_numbers = []
    report_types = []
    for rcc_file in rcc_files:
        accession_numbers.append(rcc_file.split("/")[6])
        report_types.append(rcc_file.split("/")[7])

    rcc_file_df = pd.DataFrame(
        {"accession_number": accession_numbers, "report_type": report_types}
    )
    rcc_file_unique = (
        rcc_file_df.groupby(["accession_number", "report_type"]).count().reset_index()
    )
    rcc_file_unique["has_rcc_bam_file"] = True

    self.rcc = rcc_file_unique
