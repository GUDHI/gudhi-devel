# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from gudhi.datasets import remote

import shutil
import io
import sys
import pytest

from os.path import isdir, expanduser, exists
from os import remove

def test_data_home():
    # Test get_data_home and clear_data_home on new empty folder
    empty_data_home = remote.get_data_home(data_home="empty_folder_for_test")
    assert isdir(empty_data_home)

    remote.clear_data_home(data_home=empty_data_home)
    assert not isdir(empty_data_home)

def test_fetch_remote():
    # Test fetch with a wrong checksum
    with pytest.raises(OSError):
        remote._fetch_remote("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d/spiral_2d.npy", "tmp_spiral_2d.npy", file_checksum = 'XXXXXXXXXX')
    assert not exists("tmp_spiral_2d.npy")

def _get_bunny_license_print(accept_license = False):
    capturedOutput = io.StringIO()
    # Redirect stdout
    sys.stdout = capturedOutput

    bunny_arr = remote.fetch_bunny("./tmp_for_test/bunny.npy", accept_license)
    assert bunny_arr.shape == (35947, 3)
    remove("./tmp_for_test/bunny.npy")

    # Reset redirect
    sys.stdout = sys.__stdout__
    return capturedOutput

def test_print_bunny_license():
    # Test not printing bunny.npy LICENSE when accept_license = True
    assert "" == _get_bunny_license_print(accept_license = True).getvalue()
    # Test printing bunny.LICENSE file when fetching bunny.npy with accept_license = False (default)
    with open("./tmp_for_test/bunny.LICENSE") as f:
        assert f.read().rstrip("\n") == _get_bunny_license_print().getvalue().rstrip("\n")
    shutil.rmtree("./tmp_for_test")

def test_fetch_remote_datasets_wrapped():
    # Check if gudhi_data default dir exists already
    to_be_removed = not isdir(expanduser("~/gudhi_data"))
    # Test fetch_spiral_2d and fetch_bunny wrapping functions (twice, to test case of already fetched files)
    for i in range(2):
        spiral_2d_arr = remote.fetch_spiral_2d()
        assert spiral_2d_arr.shape == (114562, 2)

        bunny_arr = remote.fetch_bunny()
        assert bunny_arr.shape == (35947, 3)

    # Check that default dir was created
    assert isdir(expanduser("~/gudhi_data"))
    # Check downloaded files
    assert exists(expanduser("~/gudhi_data/points/spiral_2d/spiral_2d.npy"))
    assert exists(expanduser("~/gudhi_data/points/bunny/bunny.npy"))
    assert exists(expanduser("~/gudhi_data/points/bunny/bunny.LICENSE"))

    # Test fetch_spiral_2d and fetch_bunny wrapping functions with data directory different from default
    spiral_2d_arr = remote.fetch_spiral_2d("./another_fetch_folder_for_test/spiral_2d.npy")
    assert spiral_2d_arr.shape == (114562, 2)

    bunny_arr = remote.fetch_bunny("./another_fetch_folder_for_test/bunny.npy")
    assert bunny_arr.shape == (35947, 3)

    assert isdir("./another_fetch_folder_for_test")
    # Check downloaded files
    assert exists("./another_fetch_folder_for_test/spiral_2d.npy")
    assert exists("./another_fetch_folder_for_test/bunny.npy")
    assert exists("./another_fetch_folder_for_test/bunny.LICENSE")

    # Remove test folders
    del spiral_2d_arr
    del bunny_arr
    if to_be_removed:
        shutil.rmtree(expanduser("~/gudhi_data"))
    shutil.rmtree("./another_fetch_folder_for_test")
