from Bio import pairwise2

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union
import pandas as pd
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import rv_histogram
import numpy as np
from memoization import cached
from tqdm import tqdm

from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
import numpy as np
from scipy.optimize import differential_evolution, Bounds
import arviz as az
from matplotlib import cm
from pymc.gp.util import plot_gp_dist
import pymc as pm

# PARSeF - Probabilistic Assignment of RNA-Seq Fragments


@cached
def local_align(ref_seq, frag):
    return pairwise2.align.localms(ref_seq, frag, 5, -4, -5, -2)


@cached
def pred_clf(clf, f1_score, f2_score, f1_gaps, f2_gaps):
    if clf is None:
        raise ValueError("Classifier is not fitted. Please fit the model first.")
    if not isinstance(clf, RandomForestClassifier):
        raise TypeError("Classifier is not a RandomForestClassifier.")
    if not hasattr(clf, "predict"):
        raise AttributeError("Classifier does not have a predict method.")
    pred = clf.predict(np.array([f1_score, f2_score, f1_gaps, f2_gaps]).reshape(1, -1))[
        0
    ]
    return pred


@dataclass
class Isoforms:
    bam_fname: Union[str, Path]
    seq_fname: Union[str, Path]
    isoform_list: list[str] = field(default_factory=list)
    transcripts: dict = field(default_factory=dict)
    sequences: dict = field(default_factory=dict)
    length_distribution: Union[rv_histogram, None] = None
    left_start_distributions: dict = field(default_factory=dict)
    right_start_distributions: dict = field(default_factory=dict)
    skip_distributions: dict = field(default_factory=dict)
    effective_lengths: dict = field(default_factory=dict)
    max_skip_prob: dict = field(default_factory=dict)
    y: Union[pd.Series, None] = None
    X: Union[pd.DataFrame, None] = None
    clf: Union[RandomForestClassifier, None] = None
    splits: list = field(default_factory=list)
    normalized_rna_seq_vector: Union[np.ndarray, None] = None
    distances: list = field(default_factory=list)
    transcript_count: int = 0
    simulation_results: Union[dict, None] = None
    simulated_data: Union[pd.DataFrame, None] = None
    sample: dict = field(default_factory=dict)
    hdi: dict = field(default_factory=dict)
    fractions: Union[np.ndarray, None] = None
    fraction_neg_cis: Union[list, None] = None
    fraction_pos_cis: Union[list, None] = None
    isoform_fractions: Union[pd.DataFrame, None] = None

    def simulate_experiment(self, n=None, max_iter=10):
        # Initial guess
        x = [1 / (i + 1) for i in range(len(self.sequences))]
        x0 = x[1:][::-1]
        print(x0)
        if n is None:
            n = self.transcript_count

        # Define bounds for each variable
        bounds = Bounds(lb=[0] * len(x0), ub=[1] * len(x0))

        # Using a sequence of tuples
        # bounds = ((-1, 2), (0, 3))
        objective_function = lambda x: self._score_x(x, n)
        # Minimize the function with bounds using L-BFGS-B method
        self.simulation_results = differential_evolution(
            objective_function, bounds, x0=x0, maxiter=max_iter
        )
        split_transcripts = list(self.sequences.keys())[:-1]
        simulated_data = {key: [] for key in split_transcripts}
        for i in range(len(self.splits)):
            for j, t in enumerate(split_transcripts):
                simulated_data[t].append(self.splits[i][j])
        simulated_data["distances"] = self.distances
        self.simulated_data = pd.DataFrame(simulated_data)

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
        save_label = ""
        with open(self.seq_fname) as fp:
            for line in fp.readlines():
                if save_next:
                    self.sequences[save_label] = line.strip()
                    save_next = False
                else:
                    for isoform in self.isoform_list:
                        if line.find(isoform) > -1:
                            save_next = True
                            save_label = isoform
                            break

    def _measure_length_distribution(self):
        lengths = []
        for k in self.transcripts.keys():
            for isoform in self.transcripts[k]:
                lengths.append(len(isoform.seq))

        length_hist = np.histogram(lengths)

        self.length_distribution = rv_histogram(length_hist, density=False)

    def _measure_transcript_distribution(self):
        self.left_start_distributions = {}
        self.right_start_distributions = {}
        self.skip_distributions = {}
        self.effective_lengths = {}
        self.max_skip_prob = {}
        for isoform_name in self.transcripts.keys():
            transcripts = self.transcripts[isoform_name]
            skips = []
            left_start = []
            right_start = []
            for transcript in transcripts:
                if transcript.is_forward:
                    left_start.append(int(transcript.to_dict()["ref_pos"]) - 1)
                    right_start.append(int(transcript.to_dict()["next_ref_pos"]) - 1)
                    skips.append(right_start[-1] - left_start[-1])

            ls_hist = np.histogram(left_start)

            ls_hist_dist = rv_histogram(ls_hist, density=False)

            rs_hist = np.histogram(right_start)
            rs_hist_dist = rv_histogram(rs_hist, density=False)

            sk_hist = np.histogram(skips)
            sk_hist_dist = rv_histogram(sk_hist, density=False)

            self.left_start_distributions[isoform_name] = ls_hist_dist
            self.right_start_distributions[isoform_name] = rs_hist_dist
            self.skip_distributions[isoform_name] = sk_hist_dist

            p = []
            for i in range(max(skips) + 1):
                p.append(sk_hist_dist.pdf(i))

            self.effective_lengths[isoform_name] = max(right_start) - min(left_start)
            self.max_skip_prob[isoform_name] = max(p)

    def _sim_int(self, hist_dist):
        return int(np.round(hist_dist.ppf(np.random.rand()) + 0.25))

    def sim_isoform_transcripts(self, isoform_counts):
        sim_left_starts = {}
        sim_right_starts = {}
        sim_skips = {}
        for transcript in isoform_counts.keys():
            sim_transcript_left_starts = []
            sim_transcript_right_starts = []
            sim_transcript_skips = []
            while len(sim_transcript_left_starts) < isoform_counts[transcript]:
                ls = self._sim_int(self.left_start_distributions[transcript])
                rs = self._sim_int(self.right_start_distributions[transcript])
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

    def sim_isoform_counts(self, n, splits: list):
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
        for_frag_len = self._sim_int(self.length_distribution)
        rev_frag_len = self._sim_int(self.length_distribution)
        return (
            seq[sim_left_start : sim_left_start + for_frag_len],
            seq[sim_right_start : sim_right_start + rev_frag_len],
        )

    def _generate_alignement_model_features(self):
        transcript_dfs = []
        all_names = []
        for key in self.transcripts.keys():
            names = []
            ori = []
            seq = []
            for read in self.transcripts[key]:
                names.append(read.to_dict()["name"])
                ori.append(read.is_forward)
                seq.append(read.seq)

            transcript_dfs.append(
                pd.DataFrame(
                    {"names": names, "ori": ori, "seq": seq, "transcript": key}
                )
            )
            all_names += names
        transcript_count = {}
        overlaps = set(all_names)
        all_names = list(overlaps)
        self.transcript_count = len(all_names)

        for transcript_df in transcript_dfs:
            transcript = transcript_df.transcript.iloc[0]
            name_set = set(transcript_df.names)
            unique_set = name_set.copy()
            for inner_transcript_df in transcript_dfs:
                if inner_transcript_df.transcript.iloc[0] == transcript:
                    continue
                inner_names = set(inner_transcript_df.names)
                unique_set = unique_set - inner_names
            transcript_count[transcript] = len(unique_set)
            overlaps = overlaps - unique_set

        for i in range(len(overlaps)):
            hits = 0
            transcript_hits = {key: 0 for key in self.transcripts.keys()}

            for transcript in self.transcripts.keys():
                for k in self.transcripts[transcript]:
                    if k.to_dict()["name"] == list(overlaps)[i]:
                        hits += 1
                        transcript_hits[transcript] += 1

            for transcript in transcript_hits.keys():
                transcript_count[transcript] += transcript_hits[transcript] / hits
        transcript_count_df = pd.Series(transcript_count)
        self.normalized_rna_seq_vector = (
            np.array(transcript_count_df.values) / transcript_count_df.sum()
        )
        dfs = []

        for name in tqdm(all_names):
            matched = []
            for_seq = ""
            rev_seq = ""
            for transcript_df in transcript_dfs:
                match_df = transcript_df[transcript_df.names == name]
                if len(match_df) == 2:
                    for_seq = match_df[match_df.ori].iloc[0].seq
                    rev_seq = match_df[match_df.ori == False].iloc[0].seq
                    matched.append(1)
                else:
                    matched.append(0)

            for_scores = []
            rev_scores = []
            for_gaps = []
            rev_gaps = []
            # for_len_mismatches = []
            # rev_len_mismatches = []
            for key in self.sequences.keys():

                l = local_align(for_seq, self.sequences[key])
                for_scores.append(l[0].score)
                for_match = l[0].seqA[l[0].start : l[0].end]
                for_gaps.append(len(for_match.split("-")) - 1)
                # for_len_mismatches.append(len(for_seq) - len(for_match))

                l = local_align(rev_seq, self.sequences[key])
                rev_scores.append(l[0].score)
                rev_match = l[0].seqA[l[0].start : l[0].end]
                rev_gaps.append(len(rev_match.split("-")) - 1)
                # rev_len_mismatches.append(len(rev_seq) - len(rev_match))

            inside_df = pd.DataFrame(
                {
                    "f_score": for_scores,
                    "r_score": rev_scores,
                    "f_gaps": for_gaps,
                    "r_gaps": rev_gaps,
                    # "f_len_mismatch": for_len_mismatches,
                    # "r_len_mismatch": rev_len_mismatches,
                    "match": matched,
                }
            )
            dfs.append(inside_df)
        matched_df = pd.concat(dfs)

        self.y = matched_df.match
        self.X = matched_df.drop("match", axis=1)

    def fit_alignment_model(self):
        if self.X is None:
            self._generate_alignement_model_features()
        assert self.X is not None, "X is None, cannot fit model."
        assert self.y is not None, "y is None, cannot fit model."
        X_train, X_test, y_train, y_test = train_test_split(
            self.X, self.y, test_size=0.33
        )
        clf = RandomForestClassifier(max_depth=10, random_state=42)
        clf.fit(X_train, y_train)

        model_accuracy = clf.score(X_test, y_test)
        print(f"Model accuracy: {model_accuracy}")
        if model_accuracy < 0.98:
            print(
                "Model accuracy is below 0.98, consider retraining with more data or adjusting parameters."
            )
            self.clf = None
        else:
            clf = RandomForestClassifier(max_depth=10, random_state=42)
            clf.fit(self.X.values, self.y)
            self.clf = clf

    def _x_to_split(self, xs):
        split = [xs[0]]
        for x in xs[1:]:
            rem = 1 - split[-1]
            split.append(rem * x + split[-1])
        return split

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

        split = self._x_to_split(x)
        print(f"Eval split {split}")
        self.splits.append(split)
        isoform_counts = self.sim_isoform_counts(n, split)

        sim_left_starts, sim_right_starts, sim_skips = self.sim_isoform_transcripts(
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

    def _gp_model(self, data, col):
        with pm.Model() as gp_model:
            p = pm.HalfCauchy("p", 3)
            n = pm.HalfCauchy("n", 3)
            M = pm.gp.mean.Constant(0)
            K = (n**2) * pm.gp.cov.ExpQuad(1, p)
            theta = pm.HalfNormal("theta", 50)
            t1_marg = pm.gp.Marginal(mean_func=M, cov_func=K)
            t1_marg.marginal_likelihood(
                "t1_likelihood",
                X=data[col].values.reshape(-1, 1),
                y=data.distances.values,
                sigma=theta,
            )

        with gp_model:
            idata = pm.sample(1000)

        az.plot_trace(idata)

        X_new = np.linspace(data[col].min(), data[col].max(), 100)

        with gp_model:
            f_pred = t1_marg.conditional("f_pred", X_new[:, None])

        with gp_model:
            gp_model_samples = pm.sample_posterior_predictive(
                idata.sel(draw=slice(0, 9)), var_names=["f_pred"]
            )

        f_pred_samples = az.extract(
            gp_model_samples, group="posterior_predictive", var_names=["f_pred"]
        )
        return X_new, f_pred_samples

    def _plot_gp_model(self, data, col, X_new, f_pred_samples):
        # plot the results
        fig = plt.figure(figsize=(12, 5))
        ax = fig.gca()

        # plot the samples from the gp posterior with samples and shading

        # f_pred_samples = az.extract(gp_model_samples, group="posterior_predictive", var_names=["f_pred"])
        plot_gp_dist(ax, samples=f_pred_samples.T, x=X_new)
        # plot the data and the true latent function
        # plt.plot(X, f_true, "dodgerblue", lw=3, label="True f")
        plt.plot(
            data[col], data.distances, "ok", ms=3, alpha=0.5, label="Observed data"
        )

        # axis labels and title
        plt.xlabel("X")
        plt.title("Posterior distribution over $f(x)$ at the observed values")
        plt.legend()

    def _sample_gp_fit(self, X_new, f_pred_samples):
        sample = []
        sample_str = 100
        for i in range(len(f_pred_samples[0].values)):
            v = []
            for j in range(len(f_pred_samples)):
                v.append(1 / f_pred_samples[j][i].values)
            df = pd.DataFrame({"x": X_new, "inv_v": v})
            df["cdf"] = df.inv_v.cumsum(axis=0) / df.inv_v.sum()
            for i in range(sample_str):
                r = np.random.rand()
                x_r = df[df.cdf > r].x.iloc[0]
                sample.append(x_r)

        plt.hist(sample)

        with pm.Model() as model:
            n = pm.Uniform("n", lower=0, upper=1)
            s = pm.Uniform("s", lower=-5, upper=5)
            x = pm.TruncatedNormal("n2", mu=n, sigma=s, lower=0, upper=1)
            l = pm.TruncatedNormal(
                "norm", mu=n, sigma=s, lower=0, upper=1, observed=sample
            )

        with model:
            sample_data = pm.sample()

        az.plot_trace(sample_data)

        hdi = az.hdi(sample_data, hdi_prob=0.95)
        return sample_data, hdi

    def load_data(self, fname):
        self.simulated_data = pd.read_csv(fname)

    def save_data(self, fname):
        if self.simulated_data is None:
            raise ValueError("No simulated data to save. Run simulation first.")
        self.simulated_data.to_csv(fname, index=False)

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
        for transcript in list(self.sequences.keys())[:-1]:
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
            # self.fraction_neg_cis.append(neg_ci)
            # self.fraction_pos_cis.append(pos_ci)

        self.isoform_fractions = pd.DataFrame(
            {
                "transcript": list(self.sequences.keys()),
                "fraction": self.fractions,
                "neg_ci": self.fraction_neg_cis,
                "pos_ci": self.fraction_pos_cis,
            }
        )
