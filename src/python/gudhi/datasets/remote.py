# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification
#   - 2022/09 Vincent Rouvreau: factorize _fetch_remote_license
#                               fetch activities dataset

from os.path import join, split, exists, expanduser
from os import makedirs, remove, environ

from urllib.request import urlretrieve
import hashlib
import shutil
from functools import lru_cache

import numpy as np

def _get_data_home(data_home=None):
    """
    Return the path of the remote datasets directory.
    This folder is used to store remotely fetched datasets.
    By default the datasets directory is set to a folder named 'gudhi_data' in the user home folder.
    Alternatively, it can be set by the 'GUDHI_DATA' environment variable.
    The '~' symbol is expanded to the user home folder.
    If the folder does not already exist, it is automatically created.

    Parameters
    ----------
    data_home : string
        The path to remote datasets directory.
        Default is `None`, meaning that the data home directory will be set to "~/gudhi_data",
        if the 'GUDHI_DATA' environment variable does not exist.

    Returns
    -------
    data_home: string
        The path to remote datasets directory.
    """
    if data_home is None:
        data_home = environ.get("GUDHI_DATA", join("~", "gudhi_data"))
    data_home = expanduser(data_home)
    makedirs(data_home, exist_ok=True)
    return data_home


def clear_data_home(data_home=None):
    """
    Delete the data home cache directory and all its content.

    Parameters
    ----------
    data_home : string, default is None.
        The path to remote datasets directory.
        If `None` and the 'GUDHI_DATA' environment variable does not exist,
        the default directory to be removed is set to "~/gudhi_data".
    """
    # On windows, needs to clear cache before removing a file - no lazy file deletion
    _load_and_cache_activity.cache_clear()
    data_home = _get_data_home(data_home)
    shutil.rmtree(data_home)


def _checksum_sha256(file_path):
    """
    Compute the file checksum using sha256.

    Parameters
    ----------
    file_path: string
        Full path of the created file including filename.

    Returns
    -------
        The hex digest of file_path.
    """
    sha256_hash = hashlib.sha256()
    chunk_size = 4096
    with open(file_path, "rb") as f:
        # Read and update hash string value in blocks of 4K
        while True:
            buffer = f.read(chunk_size)
            if not buffer:
                break
            sha256_hash.update(buffer)
    return sha256_hash.hexdigest()


def _fetch_remote(url, file_path, file_checksum=None):
    """
    Fetch the wanted dataset from the given url and save it in file_path.

    Parameters
    ----------
    url : string
        The url to fetch the dataset from.
    file_path : string
        Full path of the downloaded file including filename.
    file_checksum : string
        The file checksum using sha256 to check against the one computed on the downloaded file.
        Default is 'None', which means the checksum is not checked.

    Raises
    ------
    OSError
        If the computed SHA256 checksum of file does not match the one given by the user.
    """

    # Get the file
    urlretrieve(url, file_path)

    if file_checksum is not None:
        checksum = _checksum_sha256(file_path)
        if file_checksum != checksum:
            # Remove file and raise error
            remove(file_path)
            raise OSError(
                "{} has a SHA256 checksum : {}, "
                "different from expected : {}."
                "The file may be corrupted or the given url may be wrong !".format(file_path, checksum, file_checksum)
            )


def _fetch_remote_license(license_url, license_path, license_checksum=None, accept_license=False):
    """
    Fetch the wanted license from the given url and save it in file_path.

    Parameters
    ----------
    license_url : string
        The url to fetch the license file from.
    license_path : string
        Full path of the downloaded file including filename.
    license_checksum : string
        The license file checksum using sha256 to check against the one computed on the downloaded file.
        Default is 'None', which means the checksum is not checked.

    Raises
    ------
    OSError
        If the computed SHA256 checksum of file does not match the one given by the user.
    """
    _fetch_remote(license_url, license_path, license_checksum)
    # Print license terms unless accept_license is set to True
    if not accept_license:
        if exists(license_path):
            with open(license_path) as f:
                print(f.read())


def _get_archive_path(file_path, label):
    """
    Get archive path based on file_path given by user and label.

    Parameters
    ----------
    file_path: string
        Full path of the file to get including filename, or None.
    label: string
        Label used along with 'data_home' to get archive path, in case 'file_path' is None.

    Returns
    -------
        Full path of archive including filename.
    """
    if file_path is None:
        archive_path = join(_get_data_home(), label)
        dirname = split(archive_path)[0]
        makedirs(dirname, exist_ok=True)
    else:
        archive_path = file_path
        dirname = split(archive_path)[0]
        makedirs(dirname, exist_ok=True)

    return archive_path


def fetch_spiral_2d(file_path=None):
    """
    Load the spiral_2d dataset.

    Note that if the dataset already exists in the target location, it is not downloaded again,
    and the corresponding array is returned from cache.

    Parameters
    ----------
    file_path : string
        Full path of the downloaded file including filename.

        Default is None, meaning that it's set to "data_home/points/spiral_2d/spiral_2d.npy".

        The "data_home" directory is set by default to "~/gudhi_data",
        unless the 'GUDHI_DATA' environment variable is set.

    Returns
    -------
    points: numpy array
        Array of shape (114562, 2).
    """
    file_url = "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d/spiral_2d.npy"
    file_checksum = "2226024da76c073dd2f24b884baefbfd14928b52296df41ad2d9b9dc170f2401"

    archive_path = _get_archive_path(file_path, "points/spiral_2d/spiral_2d.npy")

    if not exists(archive_path):
        _fetch_remote(file_url, archive_path, file_checksum)

    return np.load(archive_path, mmap_mode="r")


def fetch_bunny(file_path=None, accept_license=False):
    """
    Load the Stanford bunny dataset.

    This dataset contains 35947 vertices.

    Note that if the dataset already exists in the target location, it is not downloaded again,
    and the corresponding array is returned from cache.

    Parameters
    ----------
    file_path : string
        Full path of the downloaded file including filename.

        Default is None, meaning that it's set to "data_home/points/bunny/bunny.npy".
        In this case, the LICENSE file would be downloaded as "data_home/points/bunny/bunny.LICENSE".

        The "data_home" directory is set by default to "~/gudhi_data",
        unless the 'GUDHI_DATA' environment variable is set.

    accept_license : boolean
        Flag to specify if user accepts the file LICENSE and prevents from printing the corresponding license terms.

        Default is False.

    Returns
    -------
    points: numpy array
        Array of shape (35947, 3).
    """

    file_url = "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/bunny/bunny.npy"
    file_checksum = "f382482fd89df8d6444152dc8fd454444fe597581b193fd139725a85af4a6c6e"
    license_url = "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/bunny/bunny.LICENSE"
    license_checksum = "b763dbe1b2fc6015d05cbf7bcc686412a2eb100a1f2220296e3b4a644c69633a"

    archive_path = _get_archive_path(file_path, "points/bunny/bunny.npy")

    if not exists(archive_path):
        _fetch_remote(file_url, archive_path, file_checksum)
        license_path = join(split(archive_path)[0], "bunny.LICENSE")
        _fetch_remote_license(license_url, license_path, license_checksum, accept_license)

    return np.load(archive_path, mmap_mode="r")


@lru_cache(maxsize=None)
def _load_and_cache_activity(file_path):
    return np.load(file_path, mmap_mode="r")


def fetch_daily_activities(file_path=None, subset=None, accept_license=False):
    """
    Load a subset of the Daily and Sports Activities dataset. This dataset comes from
    https://archive.ics.uci.edu/ml/datasets/daily+and+sports+activities (CC BY 4.0 license).

    Note that if the dataset already exists in the target location, it is not downloaded again,
    and the corresponding dataset is read from cache.

    Parameters
    ----------
    file_path : string
        Full path of the downloaded file including filename.

        Default is None, meaning that it's set to "data_home/points/activities/activities_p1_left_leg.npy".
        In this case, the LICENSE file would be downloaded as "data_home/points/activities/activities.LICENSE".

        The "data_home" directory is set by default to "~/gudhi_data",
        unless the 'GUDHI_DATA' environment variable is set.

    subset : string
        This argument allows to download the following subsets:
         * 'cross_training' Only left leg magnetometer of cross training activity performed by the person 1.
           It contains 7.500 vertices in dimension 3.
         * 'jumping' Only left leg magnetometer of jumping activity performed by the person 1.
           It contains 7.500 vertices in dimension 3.
         * 'stepper' Only left leg magnetometer of stepper activity performed by the person 1.
           It contains 7.500 vertices in dimension 3.
         * 'walking' Only left leg magnetometer of walking activity performed by the person 1.
           It contains 7.500 vertices in dimension 3.
         * None (default value) This dataset contains 30.000 vertices in dimension 3 + activity type column
           (`14.` for 'cross_training', `18.` for 'jumping', `13.` for 'stepper', or `9.` for 'walking')

    accept_license : boolean
        Flag to specify if user accepts the file LICENSE and prevents from printing the corresponding license terms.

        Default is False.

    Returns
    -------
    points: numpy array
        Depending on subset value:
         * Array of shape (7.500, 3, dtype = float).

        Or
         * Array of shape (30.000, 4, dtype = float) when `subset is None`.
    """
    # To avoid download when subset is not correct
    if subset not in ["walking", "stepper", "cross_training", "jumping", None]:
        raise ValueError("Unknown subset value")

    file_url = "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/activities/activities_p1_left_leg.npy"
    file_checksum = "ff813f717dbd3c8f8a95e59a7d8496d0b43c440e85a4db1f2b338cbfa9c02f25"
    activities_license_url = (
        "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/activities/activities.LICENSE"
    )
    activities_license_checksum = "f5ce6749fa9d5359d7b0c4c37d0b61e5d9520f9494cd53be94295d3967ee4023"
    gudhi_data_set_path = "points/activities/activities_p1_left_leg.npy"

    archive_path = _get_archive_path(file_path, gudhi_data_set_path)
    if not exists(archive_path):
        _fetch_remote(file_url, archive_path, file_checksum)
        license_path = join(split(archive_path)[0], "activities.LICENSE")
        _fetch_remote_license(activities_license_url, license_path, activities_license_checksum, accept_license)

    ds = _load_and_cache_activity(archive_path)

    if subset is not None:
        activity_pos = {"walking": 0, "stepper": 1, "cross_training": 2, "jumping": 3}
        per_activity = 7500
        start_pos = activity_pos[subset] * per_activity
        ds = ds[start_pos : (start_pos + per_activity), 0:3]

    return ds
