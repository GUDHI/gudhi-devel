# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification


from gudhi.datasets import remote
import re
from os.path import isfile, exists
from os import makedirs
import io
import sys
import pytest

def _check_dir_file_names(path_file_dw, filename, dirname):
    assert isfile(path_file_dw)

    names_dw = re.split(r' |/|\\', path_file_dw)
    # Case where inner directories are created in "remote_datasets/"; e.g: "remote_datasets/bunny"
    if len(names_dw) >= 3:
        for i in range(len(names_dw)-1):
            assert re.split(r' |/|\\', dirname)[i] == names_dw[i]
        assert filename == names_dw[i+1]
    else:
        assert dirname == names_dw[0]
        assert filename == names_dw[1]

def _check_fetch_output(url, filename, dirname = "remote_datasets", file_checksum = None):
    if not exists(dirname):
        makedirs(dirname)
    path_file_dw = remote._fetch_remote(url, filename, dirname, file_checksum)
    _check_dir_file_names(path_file_dw, filename, dirname)

def _get_bunny_license_print(accept_license = False):
    capturedOutput = io.StringIO()
    # Redirect stdout
    sys.stdout = capturedOutput

    if not exists("remote_datasets/bunny"):
        makedirs("remote_datasets/bunny")

    remote._fetch_remote("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/bunny/bunny.npy", "bunny.npy", "remote_datasets/bunny",
                 '13f7842ebb4b45370e50641ff28c88685703efa5faab14edf0bb7d113a965e1b', accept_license)
    # Reset redirect
    sys.stdout = sys.__stdout__
    return capturedOutput

def test_fetch_remote_datasets():
    # Test fetch with a wrong checksum
    with pytest.raises(OSError):
        _check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d/spiral_2d.npy", "spiral_2d.npy", file_checksum = 'XXXXXXXXXX')

    # Test files download from given urls with checksums provided
    _check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d/spiral_2d.npy", "spiral_2d.npy",
                                file_checksum = '88312ffd6df2e2cb2bde9c0e1f962d7d644c6f58dc369c7b377b298dacdc4eaf')

    _check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/sphere3D_pts_on_grid.off", "sphere3D_pts_on_grid.off",
                                file_checksum = '32f96d2cafb1177f0dd5e0a019b6ff5658e14a619a7815ae55ad0fc5e8bd3f88')

    # Test files download from given urls without checksums
    _check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d/spiral_2d.npy", "spiral_2d.npy")

    _check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/sphere3D_pts_on_grid.off", "sphere3D_pts_on_grid.off")

    # Test printing existing LICENSE file when fetching bunny.npy with accept_license = False (default)
    # Fetch LICENSE file
    if not exists("remote_datasets/bunny"):
        makedirs("remote_datasets/bunny")
    remote._fetch_remote("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/bunny/LICENSE", "LICENSE", "remote_datasets/bunny",
                 'b763dbe1b2fc6015d05cbf7bcc686412a2eb100a1f2220296e3b4a644c69633a')
    with open("remote_datasets/bunny/LICENSE") as f:
        assert f.read().rstrip("\n") == _get_bunny_license_print().getvalue().rstrip("\n")

    # Test not printing bunny.npy LICENSE when accept_license = True
    assert "" == _get_bunny_license_print(accept_license = True).getvalue()

    # Test fetch_spiral_2d and fetch_bunny wrapping functions (twice, to test case of already fetched files)
    for i in range(2):
        spiral_2d_arr = remote.fetch_spiral_2d()
        assert spiral_2d_arr.shape == (114562, 2)

        bunny_arr = remote.fetch_bunny()
        assert bunny_arr.shape == (35947, 3)
