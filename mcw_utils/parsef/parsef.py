from dataclasses import dataclass
from pathlib import Path
import pandas as pd
import numpy as np
from itertools import combinations

import numpy as np
from scipy.optimize import differential_evolution, Bounds

from .read_bam_and_seq import ReadBamAndSeq
from .distributions import Distribution
from .sim import SimulateIsoform
from .model_alignment import ModelAlignment
from .gp import GaussianProcess

# PARSeF - Probabilistic Assignment of RNA-Seq Fragments


@dataclass
class Parsef(
    ReadBamAndSeq, Distribution, SimulateIsoform, ModelAlignment, GaussianProcess
):
    def __post_init__(self):
        self.bam_fname = Path(self.bam_fname)
        self.seq_fname = Path(self.seq_fname)
        if not self.bam_fname.is_file():
            raise FileNotFoundError(f"BAM file {self.bam_fname} does not exist.")
        if not self.seq_fname.is_file():
            raise FileNotFoundError(f"Sequence file {self.seq_fname} does not exist.")
        self._read_isoform_bam()
        self._read_seq_str()

    def measure_distributions(self):
        self._measure_length_distribution()
        self._measure_transcript_distribution()

    def fit_alignment_model(self):
        self._generate_alignment_model_features()
        self._fit_alignment_model()

    def simulate_experiment(self, n=None, max_iter=10):
        # Initial guess
        if n is None:
            n = self.total_transcript_count
        split_transcripts = self.isoform_list[:-1]
        simulated_data = {isoform: [] for isoform in split_transcripts}
        basic_range = np.arange(0, 0.9999999, 1 / max_iter)
        distances = []
        for c in combinations(basic_range, len(self.isoform_list) - 1):
            distances.append(self._score_x(c, n))
            for i in range(len(c)):
                for j, t in enumerate(split_transcripts):
                    simulated_data[t].append(self.splits[i][j])
        simulated_data["distances"] = distances
        self.simulated_data = pd.DataFrame(simulated_data)

    def save_simulated_data(self, fname):
        if self.simulated_data is None:
            raise ValueError("No simulated data to save. Run simulation first.")
        self.simulated_data.to_csv(fname, index=False)

    def load_simulated_data(self, fname):
        self.simulated_data = pd.read_csv(fname)

    def measure_isoform_fractions(self, quantile=0.20):
        if self.simulated_data is None:
            raise ValueError("No simulated data to save. Run simulation first.")
        small_simulated_data = self.simulated_data[
            self.simulated_data.distances
            < self.simulated_data.distances.quantile(quantile)
        ]
        splits = []
        neg_cis = []
        pos_cis = []
        for transcript in self.isoform_list[:-1]:
            X_new, f_pred_samples = self._gp_model(small_simulated_data, transcript)

            self._plot_gp_model(small_simulated_data, transcript, X_new, f_pred_samples)

            sample, hdi = self._sample_gp_fit(X_new, f_pred_samples)
            splits.append(sample.posterior["n"].mean().values)
            neg_cis.append(hdi.n2.values[0])
            pos_cis.append(hdi.n2.values[1])
        neg_cis.append(neg_cis[-1])
        pos_cis.append(pos_cis[-1])

        self.fractions = np.diff([0] + splits + [1])
        splits.append(splits[-1])
        self.fraction_neg_cis = []
        self.fraction_pos_cis = []
        for split, neg_ci, pos_ci, fraction in zip(
            splits, neg_cis, pos_cis, self.fractions
        ):
            self.fraction_neg_cis.append((split - neg_ci) / split * fraction)
            self.fraction_pos_cis.append((pos_ci - split) / split * fraction)

        self.isoform_fractions = pd.DataFrame(
            {
                "transcript": self.isoform_list,
                "fraction": self.fractions,
                "neg_ci": self.fraction_neg_cis,
                "pos_ci": self.fraction_pos_cis,
            }
        )
