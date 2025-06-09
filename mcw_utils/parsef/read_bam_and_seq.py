from dataclasses import dataclass
import pysam
from .isoform import Isoform


@dataclass
class ReadBamAndSeq(Isoform):

    def _read_isoform_bam(self):
        samfile2 = pysam.AlignmentFile(str(self.bam_fname), "rb")
        self.transcripts = {}
        for isoform in self.isoform_list:
            self.transcripts[isoform] = []

        for isoform in self.isoform_list:
            for read in samfile2.fetch(isoform):
                self.transcripts[isoform].append(read)

    def _read_seq_str(self):
        save_next = False
        self.sequences = {}
        save_label = None
        with open(self.seq_fname) as fp:
            for line in fp.readlines():
                if save_next:
                    self.sequences[save_label] = line.strip()
                    save_next = False
                    save_label = None
                else:
                    for isoform in self.isoform_list:
                        if line.find(isoform) > -1:
                            save_next = True
                            save_label = isoform
                            break
