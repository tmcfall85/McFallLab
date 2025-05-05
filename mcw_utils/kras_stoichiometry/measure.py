import sys
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import pysam


def get_reverse(seq):
    rev_seq = []
    for s in seq:
        if s == "G":
            rev_seq.append("C")
        elif s == "C":
            rev_seq.append("G")
        elif s == "T":
            rev_seq.append("A")
        elif s == "A":
            rev_seq.append("T")
    return "".join(rev_seq)


@dataclass
class RNAseqRas:
    bam_dir: str
    bam_fnames: list[str] = field(default_factory=list)
    rsem_fnames: list[str] = field(default_factory=list)
    accession_numbers: list[str] = field(default_factory=list)
    unknown_kras_count: list[int] = field(default_factory=list)
    wt_kras_count: list[int] = field(default_factory=list)
    variant_kras_count: list[int] = field(default_factory=list)
    other_kras_count: list[int] = field(default_factory=list)
    nras_count: list[int] = field(default_factory=list)
    hras_count: list[int] = field(default_factory=list)
    other_codons: list[str] = field(default_factory=list)
    kras_variants: list[int] = field(default_factory=list)
    output_counts: pd.DataFrame = None

    def __post_init__(self):
        # Read all .bam files in the bam directory and append them to fnames
        self.bam_dir = Path(self.bam_dir)

        for item in self.bam_dir.iterdir():
            if item.is_file():
                if item.suffix == ".bam":
                    self.bam_fnames.append(item)
                    self.accession_numbers.append(item.name.split("_alig")[0])
                elif item.suffix == ".results":
                    self.rsem_fnames.append(item)

    def _read_kras_from_rsem(self, fname):
        df = pd.read_csv(fname, sep="\t")
        kras = "ENSG00000133703"
        kras_df = df[df.gene_id.str.startswith(kras)]
        return kras_df

    def _read_hras_from_rsem(self, fname):
        df = pd.read_csv(fname, sep="\t")
        hras = "ENSG00000276536"
        hras_df = df[df.gene_id.str.startswith(hras)]
        return hras_df

    def _read_nras_from_rsem(self, fname):
        df = pd.read_csv(fname, sep="\t")
        nras = "ENSG00000213281"
        nras_df = df[df.gene_id.str.startswith(nras)]
        return nras_df

    def measure(self):
        merged_data = pd.read_csv(self.bam_dir / "tempus_json_pdac_rcc_merged.csv")
        for fname in self.bam_fnames:
            if fname.with_suffix(".bam.bai").is_file() == False:
                # Create index file
                pysam.index(str(fname))
            accession_number = fname.name.split("_")[0]
            merged_acc_number = merged_data[merged_data.acc_num == accession_number]
            merged_acc_number_dna_report = merged_acc_number[
                merged_acc_number.report_type == "DNA"
            ]
            if len(merged_acc_number_dna_report) == 1:
                self._measure_kras_variant_stoichiometry(
                    fname, merged_acc_number_dna_report.iloc[0].kras_variants
                )
            else:
                self._measure_kras_variant_stoichiometry(fname, "")
        star_df = pd.DataFrame(
            {
                "accession_number": self.accession_numbers,
                "unknown_kras_count": self.unknown_kras_count,
                "wt_kras_count": self.wt_kras_count,
                "variant_kras_count": self.variant_kras_count,
                "other_kras_count": self.other_kras_count,
                "nras_count": self.nras_count,
                "hras_count": self.hras_count,
                "other_codons": self.other_codons,
                "kras_variants": self.kras_variants,
            }
        )

        kras_dfs = []
        accession_numbers = []
        for fname in self.rsem_fnames:
            accession_numbers.append(fname.name.split("_")[0])
            kras_dfs.append(self._read_kras_from_rsem(fname))

        kras_df = pd.concat(kras_dfs)
        kras_df["accession_number"] = accession_numbers

        hras_dfs = []
        accession_numbers = []
        for fname in self.rsem_fnames:
            accession_numbers.append(fname.name.split("_")[0])
            hras_dfs.append(self._read_hras_from_rsem(fname))

        hras_df = pd.concat(hras_dfs)
        hras_df["accession_number"] = accession_numbers

        nras_dfs = []
        accession_numbers = []
        for fname in self.rsem_fnames:
            accession_numbers.append(fname.name.split("_")[0])
            nras_dfs.append(self._read_nras_from_rsem(fname))

        nras_df = pd.concat(nras_dfs)
        nras_df["accession_number"] = accession_numbers

        ras_df = kras_df.merge(
            hras_df, on="accession_number", left_prefix="kras_", right_prefix="hras_"
        )
        ras_df = ras_df.merge(nras_df, on="accession_number", right_prefix="nras_")

        self.output_counts = star_df.merge(ras_df, on="accession_number")
        print(self.output_counts.head())

    def _measure_kras_variant_stoichiometry(self, fname: str, kras_variants: str):
        print(f"Measuring kras {kras_variants} stoichiometry in: {fname}")
        samfile = pysam.AlignmentFile(fname, "rb")
        # find valid contigs
        contig_names = set([record["SN"] for record in samfile.header["SQ"]])
        if "chr12" in contig_names:
            contig_prefix = "chr"
        else:
            contig_prefix = ""

        kras = []
        hras = []
        nras = []
        other_codons = []
        for read in samfile.fetch(f"{contig_prefix}12", 25209798, 25245384):
            kras.append(read)

        for read in samfile.fetch(f"{contig_prefix}11", 532693, 534322):
            hras.append(read)

        for read in samfile.fetch(f"{contig_prefix}1", 114708538, 114716160):
            nras.append(read)

        self.unknown_kras_count.append(len(kras))
        self.nras_count.append(len(nras))
        self.hras_count.append(len(hras))

        wt_count = 0
        variant_count = 0
        other_count = 0
        if kras_variants != "" and isinstance(kras_variants, str):
            if kras_variants.find("G12") > -1:
                gene_start = 25245349
                wt_codons = ["GGT", "GGC", "GGA", "GGG"]
                if kras_variants.find("G12D") > -1:
                    variant_codons = ["GAT", "GAC"]
                elif kras_variants.find("G12V") > -1:
                    variant_codons = ["GTT", "GTC", "GTA", "GTG"]
                elif kras_variants.find("G12R") > -1:
                    variant_codons = ["AGA", "AGG", "CGT", "CGC", "CGA", "CGG"]
                elif kras_variants.find("G12C") > -1:
                    variant_codons = ["TGT", "TGC"]
                elif kras_variants.find("G12A") > -1:
                    variant_codons = ["GCT", "GCC", "GCA", "GCG"]
                else:
                    print(f"UNIQUE VARIANT FOUND: {kras_variants}")
                    variant_codons = []
            elif kras_variants.find("G13") > -1:
                gene_start = 25245346
                wt_codons = ["GGT", "GGC", "GGA", "GGG"]
                if kras_variants.find("G12D") > -1:
                    variant_codons = ["GAT", "GAC"]
                elif kras_variants.find("G12C") > -1:
                    variant_codons = ["TGT", "TGC"]
                else:
                    print(f"UNIQUE VARIANT FOUND: {kras_variants}")
                    variant_codons = []
            elif kras_variants.find("Q61") > -1:
                gene_start = 25227341
                wt_codons = ["CAA", "CAG"]
                if kras_variants.find("Q61H") > -1:
                    variant_codons = ["CAT", "CAC"]
                elif kras_variants.find("Q61L") > -1:
                    variant_codons = ["CTT", "CTC", "CTA", "CTG"]
                elif kras_variants.find("Q61R") > -1:
                    variant_codons = ["AGA", "AGG", "CGT", "CGC", "CGA", "CGG"]
                else:
                    print(f"UNIQUE VARIANT FOUND: {kras_variants}")
                    variant_codons = []
            else:
                print(f"CODON NOT YET DEFINED: {kras_variants}")
                gene_start = 25245349
                wt_codons = []
                variant_codons = []

            for read in samfile.fetch(f"{contig_prefix}12", gene_start, gene_start + 2):
                codon = []
                for ap in read.aligned_pairs:
                    if (
                        ap[1] is not None
                        and ap[1] >= gene_start - 1
                        and ap[1] <= gene_start + 1
                    ):
                        if ap[0] is not None:
                            if read.is_forward:
                                try:
                                    codon.append(
                                        get_reverse(str(read.get_forward_sequence()))[
                                            ap[0]
                                        ]
                                    )
                                except IndexError:
                                    codon.append("?")
                            else:
                                try:
                                    codon.append(
                                        str(read.get_forward_sequence())[::-1][ap[0]]
                                    )
                                except IndexError:
                                    codon.append("?")
                        else:
                            codon.append("?")

                codon = "".join(codon)[::-1]
                if codon in wt_codons:
                    wt_count += 1
                elif codon in variant_codons:
                    variant_count += 1
                else:
                    other_count += 1
                    other_codons.append(codon)

        else:
            print(f"No reported kras variant for accession number: {fname.name}")

        self.wt_kras_count.append(wt_count)
        self.variant_kras_count.append(variant_count)
        self.other_kras_count.append(other_count)
        self.other_codons.append("_".join(other_codons))
        self.kras_variants.append(kras_variants)

    def write(self):
        self.output_counts.to_csv(self.bam_dir / "ras_output_counts.csv", index=False)


if __name__ == "__main__":
    rnaseq_ras = RNAseqRas(sys.argv[1])
    rnaseq_ras.measure()
    rnaseq_ras.write()
