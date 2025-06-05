from dataclasses import dataclass, field
from typing import Union
import pandas as pd
import numpy as np
from scipy.stats import rv_histogram
from sklearn.ensemble import RandomForestClassifier
from pathlib import Path


@dataclass
class Isoform:
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
    total_transcript_count: int = 0
    simulation_results: Union[dict, None] = None
    simulated_data: Union[pd.DataFrame, None] = None
    sample: dict = field(default_factory=dict)
    hdi: dict = field(default_factory=dict)
    fractions: Union[np.ndarray, None] = None
    fraction_neg_cis: Union[list, None] = None
    fraction_pos_cis: Union[list, None] = None
    isoform_fractions: Union[pd.DataFrame, None] = None
