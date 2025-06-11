from Bio.Align import PairwiseAligner

import numpy as np
from memoization import cached

from sklearn.ensemble import RandomForestClassifier
import numpy as np

aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 5
aligner.mismatch_score = -4
aligner.open_gap_score = -5
aligner.extend_gap_score = -2


@cached
def local_align(ref_seq, frag):
    return aligner.align(ref_seq, frag)


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


def sim_int(hist_dist):
    return int(np.round(hist_dist.ppf(np.random.rand()) + 0.25))


def x_to_split(xs):
    split = [xs[0]]
    for x in xs[1:]:
        rem = 1 - split[-1]
        split.append(rem * x + split[-1])
    return split
