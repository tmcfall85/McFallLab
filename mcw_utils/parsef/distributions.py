from dataclasses import dataclass, field

from scipy.stats import rv_histogram
import numpy as np

from .isoform import Isoform
from Bio.Align import PairwiseAligner
import math


@dataclass
class Distribution(Isoform):

    def _measure_length_distribution(self):
        lengths = []
        for isoform in self.isoform_list:
            for transcript in self.transcripts[isoform]:
                lengths.append(len(transcript.seq))

        length_hist = np.histogram(lengths)

        self.length_distribution = rv_histogram(length_hist, density=False)

    def _measure_transcript_distribution(self):
        self.left_start_distributions = {}
        self.right_start_distributions = {}
        self.skip_distributions = {}
        self.effective_lengths = {}
        self.max_skip_prob = {}
        for isoform in self.isoform_list:
            skips = []
            left_start = []
            right_start = []
            for transcript in self.transcripts[isoform]:
                if transcript.is_forward:
                    left_start.append(int(transcript.to_dict()["ref_pos"]) - 1)
                    right_start.append(int(transcript.to_dict()["next_ref_pos"]) - 1)

            self.min_sequences[isoform] = self.sequences[isoform][
                min(left_start) : max(right_start)
                + int(np.ceil(self.length_distribution.ppf(1)))
            ]
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 5
        aligner.mismatch_score = -4
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -4
        min_match_length = 25
        ref_weights = {}
        lcm = math.lcm(*[i + 1 for i in range(len(self.isoform_list))])
        for ref_key in self.isoform_list:
            ref_weight = np.ones(len(self.min_sequences[ref_key]))
            for tar_key in self.isoform_list:
                if tar_key != ref_key:
                    aligned = aligner.align(
                        self.min_sequences[ref_key], self.min_sequences[tar_key]
                    )
                    ref = aligned[0].aligned[0]
                    for ref_pair in ref:
                        ref_start, ref_finish = ref_pair
                        if ref_finish - ref_start > min_match_length:
                            ref_weight[ref_start:ref_finish] += 1
            ref_weight_dict = {}
            for i in range(len(ref_weight)):
                ref_weight_dict[i] = lcm // ref_weight[i]
            ref_weights[ref_key] = ref_weight_dict

        for isoform in self.isoform_list:
            skips = []
            left_start = []
            right_start = []
            for transcript in self.transcripts[isoform]:
                if transcript.is_forward:
                    ref_pos = int(transcript.to_dict()["ref_pos"]) - 1
                    next_ref_pos = int(transcript.to_dict()["next_ref_pos"]) - 1
                    try:
                        for i in range(int(ref_weights[isoform][ref_pos])):
                            left_start.append(ref_pos)
                        for i in range(int(ref_weights[isoform][next_ref_pos])):
                            right_start.append(next_ref_pos)
                    except:
                        print(self.effective_lengths)
                        print(isoform)
                        print(ref_pos, next_ref_pos)

                    skips.append(right_start[-1] - left_start[-1])

            if len(skips) == 0:
                print(isoform)
                raise ValueError(
                    f"No transcripts found for isoform {isoform}. "
                    "Please check the input BAM file."
                )

            ls_hist = np.histogram(left_start, bins="auto")
            ls_hist_dist = rv_histogram(ls_hist, density=False)

            rs_hist = np.histogram(right_start, bins="auto")
            rs_hist_dist = rv_histogram(rs_hist, density=False)

            sk_hist = np.histogram(skips, bins="auto")
            sk_hist_dist = rv_histogram(sk_hist, density=False)

            self.left_start_distributions[isoform] = ls_hist_dist
            self.right_start_distributions[isoform] = rs_hist_dist
            self.skip_distributions[isoform] = sk_hist_dist

            p = []
            for i in range(max(skips) + 1):
                p.append(sk_hist_dist.pdf(i))
            self.max_skip_prob[isoform] = max(p)
            self.effective_lengths[isoform] = max(right_start) - min(left_start)
