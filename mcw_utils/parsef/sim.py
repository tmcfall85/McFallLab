from dataclasses import dataclass
from .isoform import Isoform
from .helper import sim_int, x_to_split, local_align, pred_clf
import numpy as np
import pandas as pd
from tqdm import tqdm


@dataclass
class SimulateIsoform(Isoform):

    def _sim_isoform_counts(self, n, splits: list):
        isoform_counts = {isoform: 0 for isoform in self.isoform_list}

        for i in range(n):
            r = np.random.rand()
            assigned = False
            for split, key in zip(splits, self.isoform_list):
                if r <= split:
                    assigned = True
                    isoform_counts[key] += 1
                    break
            if assigned == False:
                isoform_counts[self.isoform_list[-1]] += 1
        return isoform_counts

    def _sim_isoform_transcripts(self, isoform_counts):
        sim_left_starts = {}
        sim_right_starts = {}
        sim_skips = {}
        for isoform in self.isoform_list:
            sim_transcript_left_starts = []
            sim_transcript_right_starts = []
            sim_transcript_skips = []
            while len(sim_transcript_left_starts) < isoform_counts[isoform]:
                ls = sim_int(self.left_start_distributions[isoform])
                rs = sim_int(self.right_start_distributions[isoform])
                sk = rs - ls

                if np.random.rand() * self.max_skip_prob[
                    isoform
                ] < self.skip_distributions[isoform].pdf(sk):
                    sim_transcript_left_starts.append(ls)
                    sim_transcript_right_starts.append(rs)
                    sim_transcript_skips.append(sk)
            sim_left_starts[isoform] = sim_transcript_left_starts
            sim_right_starts[isoform] = sim_transcript_right_starts
            sim_skips[isoform] = sim_transcript_skips
        return sim_left_starts, sim_right_starts, sim_skips

    def _generate_fragment_pair(self, seq, sim_left_start, sim_right_start):
        for_frag_len = sim_int(self.length_distribution)
        rev_frag_len = sim_int(self.length_distribution)
        return (
            seq[sim_left_start : sim_left_start + for_frag_len],
            seq[sim_right_start : sim_right_start + rev_frag_len],
        )

    def _score_frags(self, ref_seq, f1, f2):
        f1_score, f1_gaps = local_align(ref_seq, f1)

        f2_score, f2_gaps = local_align(ref_seq, f2)

        pred = pred_clf(self.clf, f1_score, f2_score, f1_gaps, f2_gaps)
        return pred

    def _score_x(self, split, n):

        # split = x_to_split(x)
        print(f"Eval split {split}")
        self.splits.append(split)
        isoform_counts = self._sim_isoform_counts(n, split)

        sim_left_starts, sim_right_starts, _ = self._sim_isoform_transcripts(
            isoform_counts
        )
        scores = {isoform: [] for isoform in self.isoform_list}
        for isoform in self.isoform_list:
            for sim_left_start, sim_right_start in tqdm(
                zip(sim_left_starts[isoform], sim_right_starts[isoform])
            ):
                for_frag, rev_frag = self._generate_fragment_pair(
                    self.sequences[isoform], sim_left_start, sim_right_start
                )
                for inner_isoform in self.isoform_list:
                    scores[inner_isoform].append(
                        self._score_frags(
                            self.min_sequences[inner_isoform], for_frag, rev_frag
                        )
                    )

        scores_df = pd.DataFrame(scores)

        scores_df["norm"] = scores_df.sum(axis=1)
        norm_isoforms = []
        for isoform in self.isoform_list:
            scores_df[f"{isoform}_norm"] = scores_df[isoform] / scores_df.norm
            norm_isoforms.append(f"{isoform}_norm")

        simulated_rna_seq_vector = np.array(scores_df[norm_isoforms].sum().values)
        simulated_rna_seq_vector /= simulated_rna_seq_vector.sum()

        print("experimental")
        print(self.normalized_rna_seq_vector)
        print("simulated")
        print(simulated_rna_seq_vector)
        assert isinstance(
            self.normalized_rna_seq_vector, np.ndarray
        ), "Normalized RNA seq vector is not an ndarray."
        assert len(simulated_rna_seq_vector) == len(
            self.normalized_rna_seq_vector
        ), "Simulated vector length does not match normalized RNA seq vector length."
        distance = np.linalg.norm(
            simulated_rna_seq_vector - self.normalized_rna_seq_vector
        )
        print(f"Distance: {distance}")
        # print(f"Cache: {local_align.cache_info()}")
        self.distances.append(distance)
        return distance
