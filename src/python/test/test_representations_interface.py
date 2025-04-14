""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       ???

    Copyright (C) 20?? Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "???"
__maintainer__ = ""
__copyright__ = "Copyright (C) 20?? Inria"
__license__ = "MIT"


import numpy as np

from sklearn.base import clone
from sklearn.cluster import KMeans

from gudhi.representations import (
    Atol,
    Landscape,
    Silhouette,
    BettiCurve,
    ComplexPolynomial,
    TopologicalVector,
    PersistenceImage,
    Entropy,
)

# WARNING: The following tests only check that the program runs, not what it outputs

vectorizers = {
    "atol": Atol(quantiser=KMeans(n_clusters=2, random_state=202312, n_init="auto")),
    # "betti": BettiCurve(),
}

diag1 = [
    np.array(
        [
            [0.0, np.inf],
            [0.0, 8.94427191],
            [0.0, 7.28010989],
            [0.0, 6.08276253],
            [0.0, 5.83095189],
            [0.0, 5.38516481],
            [0.0, 5.0],
        ]
    ),
    np.array([[11.0, np.inf], [6.32455532, 6.70820393]]),
    np.empty(shape=[0, 2]),
]

diag2 = [
    np.array(
        [
            [0.0, np.inf],
            [0.0, 8.94427191],
            [0.0, 7.28010989],
            [0.0, 6.08276253],
            [0.0, 5.83095189],
            [0.0, 5.38516481],
            [0.0, 5.0],
        ]
    ),
    np.array([[11.0, np.inf], [6.32455532, 6.70820393]]),
    np.array([[0.0, np.inf], [0.0, 1]]),
]

diag3 = [np.empty(shape=[0, 2])]


def test_fit():
    print(f" > Testing `fit`.")
    for name, vectorizer in vectorizers.items():
        print(f" >> Testing {name}")
        clone(vectorizer).fit(X=[diag1[0], diag2[0]])


def test_transform():
    print(f" > Testing `transform`.")
    for name, vectorizer in vectorizers.items():
        print(f" >> Testing {name}")
        clone(vectorizer).fit_transform(X=[diag1[0], diag2[0], diag3[0]])


def test_transform_empty():
    print(f" > Testing `transform_empty`.")
    for name, vectorizer in vectorizers.items():
        print(f" >> Testing {name}")
        copy_vec = clone(vectorizer).fit(X=[diag1[0], diag2[0]])
        copy_vec.transform(X=[diag3[0], diag3[0]])


def test_set_output():
    print(f" > Testing `set_output`.")
    try:
        import pandas

        for name, vectorizer in vectorizers.items():
            print(f" >> Testing {name}")
            clone(vectorizer).set_output(transform="pandas")
    except ImportError:
        print("Missing pandas, skipping set_output test")


def test_compose():
    print(f" > Testing composition with `sklearn.compose.ColumnTransformer`.")
    from sklearn.compose import ColumnTransformer

    for name, vectorizer in vectorizers.items():
        print(f" >> Testing {name}")
        ct = ColumnTransformer(
            [
                (f"{name}-0", clone(vectorizer), 0),
                (f"{name}-1", clone(vectorizer), 1),
                (f"{name}-2", clone(vectorizer), 2),
            ]
        )
        ct.fit_transform(X=[diag1, diag2])
