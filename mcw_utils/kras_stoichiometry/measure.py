import pandas as pd
import pysam
from dataclasses import dataclass, field
from pathlib import Path
import sys

def get_reverse(seq):
    rev_seq = []
    for s in seq:
        if s == 'G':
            rev_seq.append('C')
        elif s == 'C':
            rev_seq.append('G')
        elif s == 'T':
            rev_seq.append('A')
        elif s == 'A':
            rev_seq.append('T')
    return ''.join(rev_seq)


@dataclass
class RNAseqRas:
    bam_dir: str
    fnames: list[str] = field(default_factory=list)
    accession_numbers: list[str] = field(default_factory=list)
    unknown_kras_count: list[int] = field(default_factory=list)
    wt_kras_count: list[int] = field(default_factory=list)
    g12d_kras_count: list[int] = field(default_factory=list)
    other_kras_count: list[int] = field(default_factory=list)
    nras_count: list[int] = field(default_factory=list)
    hras_count: list[int] = field(default_factory=list)
    other_codons: list[str] = field(default_factory=list)
    output_counts: pd.DataFrame = None

    def __post_init__(self):
        # Read all .bam files in the bam directory and append them to fnames
        self.bam_dir = Path(self.bam_dir)
        
        for item in self.bam_dir.iterdir():
            if item.is_file():
                if item.suffix == '.bam':
                    self.fnames.append(item)
                    self.accession_numbers.append(item.name.split('_align')[0])
    

    def measure(self):
        for fname in self.fnames:
            self._measure_g12d_stoichiometry(fname)
        
        self.output_counts = pd.DataFrame({
            'accession_number': self.accession_numbers,
            'unknown_kras_count': self.unknown_kras_count,
            'wt_kras_count': self.wt_kras_count,
            'g12d_kras_count': self.g12d_kras_count,
            'other_kras_count': self.other_kras_count,
            'nras_count': self.nras_count,
            'hras_count': self.hras_count,
            'other_codons': self.other_codons
        })

    def _measure_g12d_stoichiometry(self, fname: str):
        print(f'Measuring kras stoichiometry in: {fname}')
        pysam.index(str(fname))
        samfile = pysam.AlignmentFile(fname, "rb")

        kras = []
        hras = []
        nras = []
        other_codons = []
        for read in samfile.fetch('chr12', 25209798,25245384):
            kras.append(read)

        for read in samfile.fetch('chr11', 532693, 534322):
            hras.append(read)
            
        for read in samfile.fetch('chr1', 114708538, 114716160):
            nras.append(read) 

        wt_count = 0
        g12d_count = 0
        other_count = 0

        for read in samfile.fetch('chr12', 25245349, 25245351):
            codon = []
            for ap in read.aligned_pairs:
                if ap[1] is not None and ap[1] >= 25245349-1 and ap[1] <= 25245351-1:
                    #print(ap)
                    if read.is_forward:                
                        codon.append(get_reverse(str(read.get_forward_sequence()))[ap[0]])
                    else:                
                        codon.append(str(read.get_forward_sequence())[::-1][ap[0]])
            codon = ''.join(codon)[::-1]
            if codon in ['GGT', 'GGC', 'GGA', 'GGG']:
                wt_count +=1
            elif codon in ['GAT', 'GAC']:
                g12d_count += 1
            else:
                other_count += 1
                other_codons.append(codon)
        
        self.unknown_kras_count.append(len(kras))
        self.wt_kras_count.append(wt_count)
        self.g12d_kras_count.append(g12d_count)
        self.other_kras_count.append(other_count)
        self.nras_count.append(len(nras))
        self.hras_count.append(len(hras))
        self.other_codons.append('_'.join(other_codons))
    
    def write(self):
        self.output_counts.to_csv(self.bam_dir / 'ras_output_counts.csv', index=False)


if __name__ == "__main__":
    rnaseq_ras = RNAseqRas(sys.argv[1])
    rnaseq_ras.measure()
    rnaseq_ras.write()