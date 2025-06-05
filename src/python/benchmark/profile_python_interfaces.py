import numpy as np
import timeit

from gudhi import RipsComplex, DelaunayComplex, DelaunayCechComplex, AlphaComplex
from gudhi.datasets.generators import points

try:
    from gudhi._pers_cub_low_dim import _persistence_on_a_line as persistence_on_a_line

    print("Cython/Pybind11 version")
except ImportError:
    from gudhi._pers_cub_low_dim_ext import _persistence_on_a_line as persistence_on_a_line

    print("Nanobind version")


def auto_timeit(stmt="pass", setup="pass", globals=None, min_number=1, min_time=0.2):
    timer = timeit.Timer(stmt=stmt, setup=setup, globals=globals)
    i = min_number
    while True:
        for j in 1, 2, 5:
            number = i * j
            time_taken = timer.timeit(number)
            if time_taken >= min_time:
                res = time_taken / number
                return res, number
        i *= 10


print(
    "\nAlgo; nb_points; Mean construction time (in ms); Nber of iterations for construction mean; Mean num_vertices time (in ms); Nber of iterations for num_vertices mean);"
)
for algo in [RipsComplex, DelaunayComplex, DelaunayCechComplex, AlphaComplex]:
    for nb_points in [1, 500, 1000, 5000, 10000]:  # , 50000, 100000]:
        pts = np.random.rand(nb_points, 3)
        # pts = points.torus(n_samples = nb_points, ambient_dim = 3)

        algo_str = str(algo).split("'")[1].split(".")[-1]
        result, nb_it = auto_timeit(stmt="algo(points=pts)", globals=globals(), min_number=5)
        print(f"{algo_str}; {nb_points}; {result * 1000:.4f}; {nb_it}; ", end="")

        result, nb_it = auto_timeit(
            stmt="stree.num_vertices()",
            setup="cplx = algo(points=pts); stree = cplx.create_simplex_tree()",
            globals=globals(),
            min_number=5,
        )
        print(f"{result * 1000:.4f}; {nb_it};")

print(
    "\nAlgo; nb_filt; dtype; Mean persistence time (in ms); Nber of iterations for persistence mean)"
)
for dtype in ["float32", "float64"]:
    for nb_filt in [1, 10000000]:
        filt = np.random.rand(nb_filt).astype(dtype)

        result, nb_it = auto_timeit(
            stmt="persistence_on_a_line(filt)", globals=globals(), min_number=5
        )
        print(f"persistence_on_a_line; {nb_filt}; {dtype}; {result * 1000:.4f}; {nb_it};")
