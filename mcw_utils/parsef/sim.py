from dataclasses import dataclass
from .isoform import Isoform
from .helper import sim_int, x_to_split, local_align, pred_clf
import numpy as np
import pandas as pd
from tqdm import tqdm


@dataclass
class SimulateIsoform(Isoform):

    def _sim_isoform_transcripts(self, isoform_counts):
        sim_left_starts = {}
        sim_right_starts = {}
        sim_skips = {}
        for transcript in isoform_counts.keys():
            sim_transcript_left_starts = []
            sim_transcript_right_starts = []
            sim_transcript_skips = []
            while len(sim_transcript_left_starts) < isoform_counts[transcript]:
                ls = sim_int(self.left_start_distributions[transcript])
                rs = sim_int(self.right_start_distributions[transcript])
                sk = rs - ls

                if np.random.rand() * self.max_skip_prob[
                    transcript
                ] < self.skip_distributions[transcript].pdf(sk):
                    sim_transcript_left_starts.append(ls)
                    sim_transcript_right_starts.append(rs)
                    sim_transcript_skips.append(sk)
            sim_left_starts[transcript] = sim_transcript_left_starts
            sim_right_starts[transcript] = sim_transcript_right_starts
            sim_skips[transcript] = sim_transcript_skips
        return sim_left_starts, sim_right_starts, sim_skips

    def _sim_isoform_counts(self, n, splits: list):
        keys = list(self.transcripts.keys())
        isoform_counts = {key: 0 for key in keys}

        for i in range(n):
            r = np.random.rand()
            assigned = False
            for split, key in zip(splits, keys):
                if r <= split:
                    assigned = True
                    isoform_counts[key] += 1
                    break
            if assigned == False:
                isoform_counts[keys[-1]] += 1
        return isoform_counts

    def _generate_fragment_pair(self, seq, sim_left_start, sim_right_start):
        for_frag_len = sim_int(self.length_distribution)
        rev_frag_len = sim_int(self.length_distribution)
        return (
            seq[sim_left_start : sim_left_start + for_frag_len],
            seq[sim_right_start : sim_right_start + rev_frag_len],
        )

    def _score_frags(self, ref_seq, f1, f2):
        l = local_align(ref_seq, f1)
        f1_score = l[0].score
        for_match = l[0].seqA[l[0].start : l[0].end]
        f1_gaps = len(for_match.split("-")) - 1

        l = local_align(ref_seq, f2)
        f2_score = l[0].score
        rev_match = l[0].seqA[l[0].start : l[0].end]
        f2_gaps = len(rev_match.split("-")) - 1

        pred = pred_clf(self.clf, f1_score, f2_score, f1_gaps, f2_gaps)
        return pred

    def _score_x(self, x, n):

        split = x_to_split(x)
        print(f"Eval split {split}")
        self.splits.append(split)
        isoform_counts = self._sim_isoform_counts(n, split)

        sim_left_starts, sim_right_starts, _ = self._sim_isoform_transcripts(
            isoform_counts
        )
        scores = {key: [] for key in self.sequences.keys()}
        for transcript in sim_left_starts.keys():
            for sim_left_start, sim_right_start in tqdm(
                zip(sim_left_starts[transcript], sim_right_starts[transcript])
            ):
                for_frag, rev_frag = self._generate_fragment_pair(
                    self.sequences[transcript], sim_left_start, sim_right_start
                )
                for inner_transcript in self.sequences.keys():
                    scores[inner_transcript].append(
                        self._score_frags(
                            self.sequences[inner_transcript], for_frag, rev_frag
                        )
                    )

        scores_df = pd.DataFrame(scores)

        scores_df["norm"] = scores_df.sum(axis=1)
        norm_keys = []
        for key in self.sequences.keys():
            scores_df[f"{key}_norm"] = scores_df[key] / scores_df.norm
            norm_keys.append(f"{key}_norm")

        simulated_rna_seq_vector = np.array(scores_df[norm_keys].sum().values)
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
        self.distances.append(distance)
        return distance
