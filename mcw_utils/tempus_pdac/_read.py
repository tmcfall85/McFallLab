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
    variants = []
    msi = []
    report_type = []
    bio_inf_pipeline_version = []
    tempus_id = []
    test_code = []
    specimen_count = []
    specimen_sample_category = []
    specimen_sample_site = []
    specimen_date = []
    specimen_block_id = []
    specimen_tumor_percentage = []
    signout_date = []
    tumor_mutational_burden = []
    tumor_mutation_burden_percentile = []
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

                        schema_version.append(data["metadata"]["schemaVersion"])
                        acc_num.append(file_path.parts[-1])
                        accession_id.append(data["order"]["accessionId"])
                        report_type.append(data["report"]["workflow"]["reportType"])
                        if "signout_date" in data["report"].keys():
                            signout_date.append(data["report"]["signout_date"])
                        else:
                            signout_date.append(data["report"]["signoutDate"])

                        has_kras_variant = False
                        if report_type[-1] == "DNA":
                            if data["metadata"]["schemaVersion"] in ["1.3.1", "1.3.2"]:
                                msi.append(data["results"]["msiStatus"] == "stable")
                            elif data["metadata"]["schemaVersion"] in ["v3.3.0"]:
                                msi.append(data["microsatelliteInstability"]["status"])
                            else:
                                msi.append(
                                    data["results"]["microsatelliteInstability"][
                                        "status"
                                    ]
                                )
                            if data["metadata"]["schemaVersion"] in ["v3.3.0"]:
                                tumor_mutational_burden.append(
                                    data["genomicVariants"]["tumorMutationalBurden"][
                                        "tmb"
                                    ]
                                )
                                tumor_mutation_burden_percentile.append(
                                    data["genomicVariants"]["tumorMutationalBurden"][
                                        "tmbPercentile"
                                    ]
                                )
                                spa_muts = data["genomicVariants"][
                                    "potentiallyActionable"
                                ]
                                sbrv_muts = data["genomicVariants"][
                                    "biologicallyRelevant"
                                ]
                                svus_muts = data["genomicVariants"][
                                    "unknownSignificance"
                                ]
                                mut_variants = []
                                mut_kras_variants = []

                                for variant in spa_muts["variants"]:
                                    if (
                                        variant["variantDetails"][0]["display"]
                                        == "KRAS"
                                    ):
                                        mut_kras_variants.append(
                                            variant["variantDetails"][0]["pVar"]
                                        )
                                    has_kras_variant = True
                                    mut_variants.append(
                                        f'{variant["variantDetails"][0]["display"]}:{variant["variantDetails"][0]["pVar"]}:somatic:potentially_actionable:{variant["variantDetails"][0]["allelicFraction"]}'
                                    )

                                for variant in sbrv_muts["variants"]:
                                    if (
                                        variant["variantDetails"][0]["display"]
                                        == "KRAS"
                                    ):
                                        mut_kras_variants.append(
                                            variant["variantDetails"][0]["pVar"]
                                        )
                                    has_kras_variant = True
                                    mut_variants.append(
                                        f'{variant["variantDetails"][0]["display"]}:{variant["variantDetails"][0]["pVar"]}:somatic:biologically_relevant:{variant["variantDetails"][0]["allelicFraction"]}'
                                    )

                                for variant in svus_muts["variants"]:
                                    if (
                                        variant["variantDetails"][0]["display"]
                                        == "KRAS"
                                    ):
                                        mut_kras_variants.append(
                                            variant["variantDetails"][0]["pVar"]
                                        )
                                    has_kras_variant = True
                                    mut_variants.append(
                                        f'{variant["variantDetails"][0]["display"]}:{variant["variantDetails"][0]["pVar"]}:somatic:unknown_significance:{variant["variantDetails"][0]["allelicFraction"]}'
                                    )

                            else:
                                tumor_mutational_burden.append(
                                    data["results"]["tumorMutationalBurden"]
                                )
                                tumor_mutation_burden_percentile.append(
                                    data["results"]["tumorMutationBurdenPercentile"]
                                )
                                spa_muts = data["results"][
                                    "somaticPotentiallyActionableMutations"
                                ]
                                sbrv_muts = data["results"][
                                    "somaticBiologicallyRelevantVariants"
                                ]
                                svus_muts = data["results"][
                                    "somaticVariantsOfUnknownSignificance"
                                ]
                                fv_muts = data["results"]["fusionVariants"]
                                irv_muts = data["results"]["inheritedRelevantVariants"][
                                    "values"
                                ]
                                ivus_muts = data["results"][
                                    "inheritedVariantsOfUnknownSignificance"
                                ]["values"]
                                mut_variants = []
                                mut_kras_variants = []
                                for mut in spa_muts:
                                    if mut["gene"] == "KRAS":
                                        for variant in mut["variants"]:
                                            mut_kras_variants.append(
                                                variant["mutationEffect"]
                                            )
                                        has_kras_variant = True

                                    for variant in mut["variants"]:
                                        mut_variants.append(
                                            f'{mut["gene"]}:{variant["mutationEffect"]}:somatic:potentially_actionable:{variant["allelicFraction"]}'
                                        )

                                for mut in sbrv_muts:
                                    if mut["gene"] == "KRAS":
                                        mut_kras_variants.append(mut["mutationEffect"])
                                        has_kras_variant = True
                                    mut_variants.append(
                                        f'{mut["gene"]}:{mut["mutationEffect"]}:somatic:biologically_relevant:{mut["allelicFraction"]}'
                                    )

                                for mut in svus_muts:
                                    if mut["gene"] == "KRAS":
                                        mut_kras_variants.append(mut["mutationEffect"])
                                        has_kras_variant = True
                                    mut_variants.append(
                                        f'{mut["gene"]}:{mut["mutationEffect"]}:somatic:unknown_significance:{mut["allelicFraction"]}'
                                    )

                                for mut in fv_muts:
                                    mut_variants.append(
                                        f'gene5={mut["gene5"]}-gene3={mut["gene3"]}:{mut["variantDescription"]}:fusion:unknown_significance:'
                                    )

                                for mut in irv_muts:
                                    if mut["gene"] == "KRAS":
                                        mut_kras_variants.append(mut["mutationEffect"])
                                        has_kras_variant = True
                                    if "allelicFraction" in mut.keys():
                                        allelic_fraction = mut["allelicFraction"]
                                    else:
                                        allelic_fraction = ""
                                    mut_variants.append(
                                        f'{mut["gene"]}:{mut["mutationEffect"]}:inherited:biologically_relevant:{allelic_fraction}'
                                    )

                                for mut in ivus_muts:
                                    if mut["gene"] == "KRAS":
                                        mut_kras_variants.append(mut["mutationEffect"])
                                        has_kras_variant = True
                                    mut_variants.append(
                                        f'{mut["gene"]}:{mut["mutationEffect"]}:inherited:unknown_significance:{mut["allelicFraction"]}'
                                    )
                            if len(mut_variants) == 0:
                                mut_variants.append("none:none:none:none:none")
                            variants.append("|".join(mut_variants))

                        else:
                            msi.append(None)
                            variants.append("")
                            tumor_mutational_burden.append(None)
                            tumor_mutation_burden_percentile.append(None)

                        if has_kras_variant == False:
                            kras_variants.append("")
                        else:
                            kras_variants.append("|".join(mut_kras_variants))

                        if data["metadata"]["schemaVersion"] in ["v3.3.0"]:
                            bio_inf_pipeline_version.append(
                                data["report"]["bioInfoPipeline"]
                            )

                        else:
                            bio_inf_pipeline_version.append(
                                data["report"]["bioInfPipeline"]
                            )

                        report_id.append(data["report"]["reportId"])
                        tempus_id.append(data["patient"]["tempusId"])
                        test_code.append(data["order"]["test"]["code"])
                        specimen_count.append(len(data["specimens"]))
                        has_tumor_specimen = False
                        for specimen in data["specimens"]:
                            # one patient has two tumor specimens which is why i put this "and" logic in
                            # i dont know this code will just take the first sample for that person
                            if (
                                specimen["sampleCategory"] == "tumor"
                                and has_tumor_specimen == False
                            ):

                                has_tumor_specimen = True
                                specimen_sample_category.append(
                                    specimen["sampleCategory"]
                                )
                                if data["metadata"]["schemaVersion"] in ["v3.3.0"]:
                                    specimen_sample_site.append(
                                        specimen["primarySampleSite"]
                                    )

                                else:
                                    specimen_sample_site.append(specimen["sampleSite"])
                                specimen_date.append(specimen["collectionDate"])
                                specimen_block_id.append(
                                    specimen["institutionData"]["blockId"]
                                )
                                specimen_tumor_percentage.append(
                                    specimen["institutionData"]["tumorPercentage"]
                                )
                        if has_tumor_specimen == False:
                            specimen_sample_category.append(
                                data["specimens"][0]["sampleCategory"]
                            )
                            if data["metadata"]["schemaVersion"] in ["v3.3.0"]:
                                specimen_sample_site.append(
                                    data["specimens"][0]["primarySampleSite"]
                                )

                            else:
                                specimen_sample_site.append(
                                    data["specimens"][0]["sampleSite"]
                                )
                            specimen_date.append(data["specimens"][0]["collectionDate"])
                            specimen_block_id.append(
                                data["specimens"][0]["institutionData"]["blockId"]
                            )
                            specimen_tumor_percentage.append(
                                data["specimens"][0]["institutionData"][
                                    "tumorPercentage"
                                ]
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
            "variants": variants,
            "msi": msi,
            "specimen_count": specimen_count,
            "specimen_sample_category": specimen_sample_category,
            "specimen_sample_site": specimen_sample_site,
            "specimen_date": specimen_date,
            "specimen_block_id": specimen_block_id,
            "specimen_tumor_percentage": specimen_tumor_percentage,
            "signout_date": signout_date,
            "tumor_mutational_burden": tumor_mutational_burden,
            "tumor_mutation_burden_percentile": tumor_mutation_burden_percentile,
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
