from dataclasses import dataclass, field

from scipy.stats import rv_histogram
import numpy as np

from .isoform import Isoform


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
                    skips.append(right_start[-1] - left_start[-1])
            if len(skips) == 0:
                print(isoform)
                raise ValueError(
                    f"No transcripts found for isoform {isoform}. "
                    "Please check the input BAM file."
                )

            ls_hist = np.histogram(
                left_start,
                bins=self.sequences[isoform] // 25,
                range=(0, len(self.sequences[isoform])),
            )

            ls_hist_dist = rv_histogram(ls_hist, density=False)

            rs_hist = np.histogram(
                right_start,
                bins=self.sequences[isoform] // 25,
                range=(0, len(self.sequences[isoform])),
            )
            rs_hist_dist = rv_histogram(rs_hist, density=False)

            sk_hist = np.histogram(skips)
            sk_hist_dist = rv_histogram(sk_hist, density=False)

            self.left_start_distributions[isoform] = ls_hist
            self.right_start_distributions[isoform] = rs_hist
            self.skip_distributions[isoform] = sk_hist_dist

            p = []
            for i in range(max(skips) + 1):
                p.append(sk_hist_dist.pdf(i))

            self.effective_lengths[isoform] = max(right_start) - min(left_start)
            self.min_sequences[isoform] = self.sequences[isoform][
                min(left_start) : max(right_start)
                + int(np.ceil(self.length_distribution.ppf(1)))
            ]
            self.max_skip_prob[isoform] = max(p)
