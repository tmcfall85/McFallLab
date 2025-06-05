from Bio import pairwise2

import numpy as np
from memoization import cached

from sklearn.ensemble import RandomForestClassifier
import numpy as np


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


def sim_int(hist_dist):
    return int(np.round(hist_dist.ppf(np.random.rand()) + 0.25))

def x_to_split(xs):
    split = [xs[0]]
    for x in xs[1:]:
        rem = 1 - split[-1]
        split.append(rem * x + split[-1])
    return split
