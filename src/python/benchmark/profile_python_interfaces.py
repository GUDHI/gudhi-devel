from gudhi import RipsComplex, DelaunayComplex, DelaunayCechComplex, AlphaComplex
from gudhi.datasets.generators import points
try:
    from gudhi._pers_cub_low_dim import _persistence_on_a_line as persistence_on_a_line
    print("Cython/Pybind11 version")
except ImportError:
    from gudhi._pers_cub_low_dim_ext import _persistence_on_a_line as persistence_on_a_line
    print("Nanobind version")

import numpy as np

import timeit

print("\nAlgo; nb_points; Construction time (in ms); num_vertices time (in ms)")
for algo in [RipsComplex, DelaunayComplex, DelaunayCechComplex, AlphaComplex]:
    for nb_points in [ 1, 500, 1000, 5000, 10000]: #, 50000, 100000]:
        pts = np.random.rand(nb_points, 3)
        # pts = points.torus(n_samples = nb_points, ambient_dim = 3)

        algo_str = str(algo).split("'")[1].split(".")[-1]
        result = timeit.timeit('algo(points=pts)', globals=globals(), number=5) * 1000
        print(f"{algo_str}; {nb_points}; {result:.4f};", end="")

        result = timeit.timeit('stree.num_vertices()',
                               setup='cplx = algo(points=pts); stree = cplx.create_simplex_tree()',
                               globals=globals(), number=5) * 1000
        print(f"{result:.4f};")

print("\nAlgo; nb_filt; dtype; Persistence time (in ms)")
for dtype in ["float32", "float64"]:
    for nb_filt in [ 1, 10000000]:
        filt = np.random.rand(nb_filt).astype(dtype)

        result = timeit.timeit('persistence_on_a_line(filt)', globals=globals(), number=5) * 1000
        print(f"persistence_on_a_line; {nb_filt}; {dtype}; {result:.4f};")
