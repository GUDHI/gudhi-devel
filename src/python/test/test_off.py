""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Marc Glisse

    Copyright (C) 2022 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Marc Glisse"
__maintainer__ = ""
__copyright__ = "Copyright (C) 2022 Inria"
__license__ = "MIT"


import numpy as np
import pytest
from tempfile import NamedTemporaryFile

import gudhi as gd


def test_off_rw():
    for dim in range(2, 6):
        X = np.random.rand(123, dim)
        gd.write_points_to_off_file("rand.off", X)
        Y = gd.read_points_from_off_file("rand.off")
        assert Y == pytest.approx(X)


def test_human_off():
    pts = gd.read_points_from_off_file("human.off")
    # Should not try to read faces
    assert pts.shape == (4706, 3)


def test_invalid_off_file():
    name = NamedTemporaryFile().name
    with open(name, "w") as f:
        f.write("nOFL\n1000 200 2\n1.2 1.3 1.4")
    with pytest.raises(ValueError):
        gd.read_points_from_off_file(name)
    # Try to open a non-existing file - a new temp file name should not exist
    with pytest.raises(FileNotFoundError):
        gd.read_points_from_off_file(NamedTemporaryFile().name)


def test_off_file_header():
    name = NamedTemporaryFile().name
    print(name)
    # OFF
    ## classic case
    random_nb_points = np.random.randint(0, 100)
    with open(name, "w") as f:
        f.write(f"OFF\n{random_nb_points} 200 2\n1.2 1.3 1.4")
    with open(name) as f:
        assert gd.off_utils._read_off_file_header(f) == (3, random_nb_points)
    ## with comments
    random_nb_points = np.random.randint(0, 100)
    with open(name, "w") as f:
        f.write(
            f"# comment on the first line\nOFF\n# comment on the third line\n{random_nb_points} 200 2\n"
            "# comment before points\n1.2 1.3 1.4"
        )
    with open(name) as f:
        assert gd.off_utils._read_off_file_header(f) == (3, random_nb_points)
    # nOFF
    ## when 'dim nb_vertices nb_faces nb_edges' on the same line
    random_nb_points = np.random.randint(0, 100)
    random_dim = np.random.randint(3, 100)
    with open(name, "w") as f:
        f.write(f"nOFF\n{random_dim} {random_nb_points} 200 2\n1.2 1.3 1.4")
    with open(name) as f:
        assert gd.off_utils._read_off_file_header(f) == (random_dim, random_nb_points)
    ## when 'dim nb_vertices nb_faces nb_edges' on the same line + comments
    random_nb_points = np.random.randint(0, 100)
    random_dim = np.random.randint(3, 100)
    with open(name, "w") as f:
        f.write(
            f"# comment on the first line\nnOFF\n# comment on the third line\n{random_dim} {random_nb_points} 200 2\n"
            "# comment before points\n1.2 1.3 1.4"
        )
    with open(name) as f:
        assert gd.off_utils._read_off_file_header(f) == (random_dim, random_nb_points)
    ## when 'dim' and 'nb_vertices nb_faces nb_edges' on the separated lines
    random_nb_points = np.random.randint(0, 100)
    random_dim = np.random.randint(3, 100)
    with open(name, "w") as f:
        f.write(f"nOFF\n{random_dim}\n{random_nb_points} 200 2\n1.2 1.3 1.4")
    with open(name) as f:
        assert gd.off_utils._read_off_file_header(f) == (random_dim, random_nb_points)
    ## when 'dim' and 'nb_vertices nb_faces nb_edges' on the separated lines + comments
    random_nb_points = np.random.randint(0, 100)
    random_dim = np.random.randint(3, 100)
    with open(name, "w") as f:
        f.write(
            f"# first comment\nnOFF\n# second comment\n{random_dim}\n# third comment\n{random_nb_points} 200 2\n"
            "# another comment\n1.2 1.3 1.4"
        )
    with open(name) as f:
        assert gd.off_utils._read_off_file_header(f) == (random_dim, random_nb_points)


def test_get_next_line():
    for end_of_line in ["\r\n", "\n"]:
        for comment_char in ["#", "@"]:
            name = NamedTemporaryFile().name
            print(name)
            with open(name, "w") as f:
                f.write(
                    f"{comment_char}{end_of_line}{end_of_line}{comment_char}{end_of_line}First{end_of_line}"
                )
                f.write(
                    f"{end_of_line}{comment_char}{end_of_line}{end_of_line}Second{end_of_line}{end_of_line}"
                )
            with open(name) as f:
                assert gd.off_utils._get_next_line(f, comment=comment_char).split() == [
                    "First"
                ]
                assert gd.off_utils._get_next_line(f, comment=comment_char).split() == [
                    "Second"
                ]
