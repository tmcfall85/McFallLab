from dataclasses import dataclass

from .isoform import Isoform
import pandas as pd
import numpy as np
from tqdm import tqdm

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import numpy as np

from .helper import local_align


@dataclass
class ModelAlignment(Isoform):

    def _generate_alignment_model_features(self):
        transcript_dfs = {}
        all_names = []
        for isoform in self.isoform_list:
            names = []
            ori = []
            seq = []
            for read in self.transcripts[isoform]:
                names.append(read.to_dict()["name"])
                ori.append(read.is_forward)
                seq.append(read.seq)

            transcript_dfs[isoform] = pd.DataFrame(
                {"names": names, "ori": ori, "seq": seq}
            )
            all_names += names
        transcript_count = {}
        overlaps = set(all_names)
        all_names = list(set(all_names))
        self.total_transcript_count = len(all_names)

        for isoform in self.isoform_list:
            unique_set = set(transcript_dfs[isoform].names)
            for inner_isoform in self.isoform_list:
                if inner_isoform == isoform:
                    continue
                inner_transcript_df = transcript_dfs[inner_isoform]
                inner_name_set = set(inner_transcript_df.names)
                unique_set = unique_set - inner_name_set
            transcript_count[isoform] = len(unique_set)
            overlaps = overlaps - unique_set

        overlaps = list(overlaps)
        for overlap in overlaps:
            hits = 0
            transcript_hits = {isoform: 0 for isoform in self.isoform_list}

            for isoform in self.isoform_list:
                for transcript in self.transcripts[isoform]:
                    if transcript.to_dict()["name"] == overlap:
                        hits += 1
                        transcript_hits[isoform] += 1

            for isoform in self.isoform_list:
                transcript_count[isoform] += transcript_hits[isoform] / hits
        transcript_count_df = pd.Series(transcript_count)
        self.normalized_rna_seq_vector = (
            np.array(transcript_count_df.values) / transcript_count_df.sum()
        )
        dfs = []

        for name in tqdm(all_names):
            matched = []
            for_seq = ""
            rev_seq = ""
            for isoform in self.isoform_list:
                transcript_df = transcript_dfs[isoform]
                match_df = transcript_df[transcript_df.names == name]
                if len(match_df) >= 2:
                    for_seq = match_df[match_df.ori].iloc[0].seq
                    rev_seq = match_df[match_df.ori == False].iloc[0].seq
                    matched.append(1)
                else:
                    matched.append(0)

            for_scores = []
            rev_scores = []
            for_gaps = []
            rev_gaps = []
            for isoform in self.isoform_list:
                for_score, for_gap = local_align(for_seq, self.min_sequences[isoform])
                rev_score, rev_gap = local_align(rev_seq, self.min_sequences[isoform])
                for_scores.append(for_score)
                rev_scores.append(rev_score)
                for_gaps.append(for_gap)
                rev_gaps.append(rev_gap)

            inside_df = pd.DataFrame(
                {
                    "f_score": for_scores,
                    "r_score": rev_scores,
                    "f_gaps": for_gaps,
                    "r_gaps": rev_gaps,
                    "match": matched,
                }
            )
            dfs.append(inside_df)
        matched_df = pd.concat(dfs)

        self.y = matched_df.match
        self.X = matched_df.drop("match", axis=1)

    def _fit_alignment_model(self):
        assert self.X is not None, "X is None, cannot fit model."
        assert self.y is not None, "y is None, cannot fit model."
        X_train, X_test, y_train, y_test = train_test_split(
            self.X, self.y, test_size=0.33
        )
        clf = RandomForestClassifier(max_depth=10, random_state=42)
        clf.fit(X_train, y_train)

        model_accuracy = clf.score(X_test, y_test)
        print(f"Model accuracy: {model_accuracy}")
        clf = RandomForestClassifier(max_depth=10, random_state=42)
        clf.fit(self.X.values, self.y)
        self.clf = clf
