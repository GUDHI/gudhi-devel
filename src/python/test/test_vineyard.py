""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Hannah Schreiber

    Copyright (C) 2026 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from pathlib import Path
import glob
import numpy as np

from gudhi import RipsComplex
from gudhi.vineyard import Vineyard, PointCloudRipsVineyard


# not really safe, but as giving command line arguments to pytest while using ctest ask for some
# more complicated gymnastic, I opted for this option for now...
def get_file_path_prefix():
    parent_dir_name = Path(__file__).resolve().parent.parent.parent.name
    if parent_dir_name == "src":
        # devel version
        return "../../../data/points/unknotting_rings/rings"
    # user version
    return "../../data/points/unknotting_rings/rings"


def test_rips_vineyard():
    path_prefix = get_file_path_prefix()
    path_suffix = ".txt"

    number_of_updates = len(glob.glob(path_prefix + "*" + path_suffix)) - 1

    rvy = PointCloudRipsVineyard.from_files(path_prefix)
    vineyard = rvy.get_current_vineyard_view()

    assert len(vineyard) != 0
    assert vineyard[0].shape[1] == number_of_updates + 1

    for step in range(0, number_of_updates + 1):
        points = np.loadtxt(path_prefix + str(step) + path_suffix)
        st = RipsComplex(points=points).create_simplex_tree(max_dimension=2)

        pers = st.persistence(homology_coeff_field=2, min_persistence=-1)
        pers = [bar[1] for bar in pers if bar[0] == 1]
        assert len(vineyard) != 1 or len(pers) == 0

        v_pers = [(b, d) for b, d in vineyard[1][:, step, :]]
        assert len(pers) == len(v_pers)

        pers.sort()
        v_pers.sort()
        assert pers == v_pers
