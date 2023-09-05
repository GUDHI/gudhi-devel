# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification
#   - 2022/09 Vincent Rouvreau: Factorize _capture_license_output to be usable by activities
#                               Test fetch activities (except pandas file) and its license file

from gudhi.datasets import remote

import shutil
import io
import sys
import pytest
import numpy as np

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
        remote._fetch_remote(
            "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d/spiral_2d.npy",
            "tmp_spiral_2d.npy",
            file_checksum="XXXXXXXXXX",
        )
    assert not exists("tmp_spiral_2d.npy")


def _capture_license_output(fetch_method, **kwargs):
    capturedOutput = io.StringIO()
    # Redirect stdout
    sys.stdout = capturedOutput

    # Force file_path in forwarded arguments
    kwargs["file_path"] = "./tmp_for_test/data.npy"
    data_array = fetch_method(**kwargs)
    # No need for basic usage, but necessary because the file is removed in the test
    del data_array
    remote._load_and_cache_activity.cache_clear()
    # No need to keep numpy file
    remove("./tmp_for_test/data.npy")

    # Reset redirect
    sys.stdout = sys.__stdout__
    return capturedOutput


def test_bunny_license_capture():
    # Test not printing bunny.npy LICENSE when accept_license = True
    assert "" == _capture_license_output(remote.fetch_bunny, accept_license=True).getvalue()
    # Test printing bunny.LICENSE file when fetching bunny.npy with accept_license = False (default)
    with open("./tmp_for_test/bunny.LICENSE") as f:
        assert f.read().rstrip("\n") == _capture_license_output(remote.fetch_bunny).getvalue().rstrip("\n")
    shutil.rmtree("./tmp_for_test")


def test_activities_license_capture():
    for subset in [None, "cross_training", "jumping", "stepper", "walking"]:
        # Test not printing activities.LICENSE when accept_license = True
        assert (
            "" == _capture_license_output(remote.fetch_daily_activities, subset=subset, accept_license=True).getvalue()
        )
        # Test printing activities.LICENSE file when fetching "activities".npy with accept_license = False (default)
        with open("./tmp_for_test/activities.LICENSE") as f:
            assert f.read().rstrip("\n") == _capture_license_output(
                remote.fetch_daily_activities, subset=subset
            ).getvalue().rstrip("\n")
        shutil.rmtree("./tmp_for_test")


def test_activities_unknown_subset():
    with pytest.raises(ValueError):
        remote.fetch_daily_activities(subset="some weird subset value")


def test_fetch_remote_datasets_wrapped():
    # Test fetch_spiral_2d and fetch_bunny wrapping functions with data directory different from default (twice, to test case of already fetched files)
    # Default case is not tested because it would fail in case the user sets the 'GUDHI_DATA' environment variable locally
    for i in range(2):
        spiral_2d_arr = remote.fetch_spiral_2d("./another_fetch_folder_for_test/spiral_2d.npy")
        assert spiral_2d_arr.shape == (114562, 2)
        assert spiral_2d_arr.dtype == np.dtype('float32')

        bunny_arr = remote.fetch_bunny("./another_fetch_folder_for_test/bunny.npy")
        assert bunny_arr.shape == (35947, 3)
        assert bunny_arr.dtype == np.dtype('float32')

        activities_arr = remote.fetch_daily_activities("./another_fetch_folder_for_test/activities.npy", subset=None)
        assert activities_arr.shape == (30000, 4)
        assert activities_arr.dtype == np.dtype('float32')

        # Order is important
        walking_arr = remote.fetch_daily_activities("./another_fetch_folder_for_test/activities.npy", subset="walking")
        np.testing.assert_array_equal(walking_arr, activities_arr[activities_arr[:,3] == 9.][:,0:3])
        stepper_arr = remote.fetch_daily_activities("./another_fetch_folder_for_test/activities.npy", subset="stepper")
        np.testing.assert_array_equal(stepper_arr, activities_arr[activities_arr[:,3] == 13.][:,0:3])
        cross_training_arr = remote.fetch_daily_activities("./another_fetch_folder_for_test/activities.npy", subset="cross_training")
        np.testing.assert_array_equal(cross_training_arr, activities_arr[activities_arr[:,3] == 14.][:,0:3])
        jumping_arr = remote.fetch_daily_activities("./another_fetch_folder_for_test/activities.npy", subset="jumping")
        np.testing.assert_array_equal(jumping_arr, activities_arr[activities_arr[:,3] == 18.][:,0:3])
    
    # No need for basic usage, but necessary because the file is removed in the test
    del spiral_2d_arr, bunny_arr, activities_arr, walking_arr, stepper_arr, cross_training_arr, jumping_arr
    remote._load_and_cache_activity.cache_clear()

    # Check that the directory was created
    assert isdir("./another_fetch_folder_for_test")
    # Check downloaded files
    assert exists("./another_fetch_folder_for_test/spiral_2d.npy")
    assert exists("./another_fetch_folder_for_test/bunny.npy")
    assert exists("./another_fetch_folder_for_test/bunny.LICENSE")

    assert exists("./another_fetch_folder_for_test/activities.npy")
    assert exists("./another_fetch_folder_for_test/activities.LICENSE")

    # Remove test folders
    shutil.rmtree("./another_fetch_folder_for_test")


def test_gudhi_data_env():
    # Set environment variable "GUDHI_DATA"
    environ["GUDHI_DATA"] = "./test_folder_from_env_var"
    bunny_arr = remote.fetch_bunny()
    assert exists("./test_folder_from_env_var/points/bunny/bunny.npy")
    assert exists("./test_folder_from_env_var/points/bunny/bunny.LICENSE")
    # No need for basic usage, but necessary because the file is removed in the test
    del bunny_arr
    
    # Remove test folder
    shutil.rmtree("./test_folder_from_env_var")
