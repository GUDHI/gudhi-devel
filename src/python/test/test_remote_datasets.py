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
import shutil
import io
import sys
import pytest

from os.path import isfile, isdir, expanduser
from os import makedirs

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
    makedirs(dirname, exist_ok=True)
    path_file_dw = remote._fetch_remote(url, filename, dirname, file_checksum)
    _check_dir_file_names(path_file_dw, filename, dirname)

def _get_bunny_license_print(accept_license = False):
    capturedOutput = io.StringIO()
    # Redirect stdout
    sys.stdout = capturedOutput

    makedirs("remote_datasets/bunny", exist_ok=True)

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
    makedirs("remote_datasets/bunny", exist_ok=True)
    remote._fetch_remote("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/bunny/LICENSE", "LICENSE", "remote_datasets/bunny",
                 'b763dbe1b2fc6015d05cbf7bcc686412a2eb100a1f2220296e3b4a644c69633a')
    with open("remote_datasets/bunny/LICENSE") as f:
        assert f.read().rstrip("\n") == _get_bunny_license_print().getvalue().rstrip("\n")

    # Test not printing bunny.npy LICENSE when accept_license = True
    assert "" == _get_bunny_license_print(accept_license = True).getvalue()

    # Remove "remote_datasets" directory and all its content
    shutil.rmtree("remote_datasets")

def test_fetch_remote_datasets_wrapped():
    # Test fetch_spiral_2d and fetch_bunny wrapping functions (twice, to test case of already fetched files)
    for i in range(2):
        spiral_2d_arr = remote.fetch_spiral_2d()
        assert spiral_2d_arr.shape == (114562, 2)

        bunny_arr = remote.fetch_bunny()
        assert bunny_arr.shape == (35947, 3)

    # Check that default dir was created
    assert isdir(expanduser("~/remote_datasets"))

    # Test fetch_spiral_2d and fetch_bunny wrapping functions with data directory different from default
    spiral_2d_arr = remote.fetch_spiral_2d(dirname = "~/another_fetch_folder")
    assert spiral_2d_arr.shape == (114562, 2)

    bunny_arr = remote.fetch_bunny(dirname = "~/another_fetch_folder")
    assert bunny_arr.shape == (35947, 3)

    assert isdir(expanduser("~/another_fetch_folder"))

    # Remove test folders
    del spiral_2d_arr
    del bunny_arr
    shutil.rmtree(expanduser("~/remote_datasets"))
    shutil.rmtree(expanduser("~/another_fetch_folder"))

    assert not isdir(expanduser("~/remote_datasets"))
    assert not isdir(expanduser("~/another_fetch_folder"))

def test_data_home():
    # Test get_data_home and clear_data_home on new empty folder
    empty_data_home = remote.get_data_home(data_home="empty_folder")
    assert isdir(empty_data_home)

    remote.clear_data_home(data_home=empty_data_home)
    assert not isdir(empty_data_home)
