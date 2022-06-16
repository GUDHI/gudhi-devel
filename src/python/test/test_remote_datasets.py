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
from os import remove, environ

def test_data_home():
    # Test _get_data_home and clear_data_home on new empty folder
    empty_data_home = remote._get_data_home(data_home="empty_folder_for_test")
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
    del bunny_arr
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
    # Test fetch_spiral_2d and fetch_bunny wrapping functions with data directory different from default (twice, to test case of already fetched files)
    # Default case is not tested because it would fail in case the user sets the 'GUDHI_DATA' environment variable locally
    for i in range(2):
        spiral_2d_arr = remote.fetch_spiral_2d("./another_fetch_folder_for_test/spiral_2d.npy")
        assert spiral_2d_arr.shape == (114562, 2)

        bunny_arr = remote.fetch_bunny("./another_fetch_folder_for_test/bunny.npy")
        assert bunny_arr.shape == (35947, 3)

    # Check that the directory was created
    assert isdir("./another_fetch_folder_for_test")
    # Check downloaded files
    assert exists("./another_fetch_folder_for_test/spiral_2d.npy")
    assert exists("./another_fetch_folder_for_test/bunny.npy")
    assert exists("./another_fetch_folder_for_test/bunny.LICENSE")

    # Remove test folders
    del spiral_2d_arr
    del bunny_arr
    shutil.rmtree("./another_fetch_folder_for_test")

def test_gudhi_data_env():
    # Set environment variable "GUDHI_DATA"
    environ["GUDHI_DATA"] = "./test_folder_from_env_var"
    bunny_arr = remote.fetch_bunny()
    assert bunny_arr.shape == (35947, 3)
    assert exists("./test_folder_from_env_var/points/bunny/bunny.npy")
    assert exists("./test_folder_from_env_var/points/bunny/bunny.LICENSE")
    # Remove test folder
    del bunny_arr
    shutil.rmtree("./test_folder_from_env_var")
