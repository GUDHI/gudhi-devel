# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from os.path import join, split, exists, expanduser
from os import makedirs, remove, environ

from urllib.request import urlretrieve
import hashlib
import shutil

import numpy as np

def _get_data_home(data_home = None):
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


def clear_data_home(data_home = None):
    """
    Delete the data home cache directory and all its content.

    Parameters
    ----------
    data_home : string, default is None.
        The path to remote datasets directory.
        If `None` and the 'GUDHI_DATA' environment variable does not exist,
        the default directory to be removed is set to "~/gudhi_data".
    """
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
    with open(file_path,"rb") as f:
        # Read and update hash string value in blocks of 4K
        while True:
            buffer = f.read(chunk_size)
            if not buffer:
                break
            sha256_hash.update(buffer)
    return sha256_hash.hexdigest()

def _fetch_remote(url, file_path, file_checksum = None):
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
    IOError
        If the computed SHA256 checksum of file does not match the one given by the user.
    """

    # Get the file
    urlretrieve(url, file_path)

    if file_checksum is not None:
        checksum = _checksum_sha256(file_path)
        if file_checksum != checksum:
            # Remove file and raise error
            remove(file_path)
            raise IOError("{} has a SHA256 checksum : {}, "
                        "different from expected : {}."
                        "The file may be corrupted or the given url may be wrong !".format(file_path, checksum, file_checksum))

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

def fetch_spiral_2d(file_path = None):
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
    file_checksum = '2226024da76c073dd2f24b884baefbfd14928b52296df41ad2d9b9dc170f2401'

    archive_path = _get_archive_path(file_path, "points/spiral_2d/spiral_2d.npy")

    if not exists(archive_path):
        _fetch_remote(file_url, archive_path, file_checksum)

    return np.load(archive_path, mmap_mode='r')

def fetch_bunny(file_path = None, accept_license = False):
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
    file_checksum = 'f382482fd89df8d6444152dc8fd454444fe597581b193fd139725a85af4a6c6e'
    license_url = "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/bunny/bunny.LICENSE"
    license_checksum = 'b763dbe1b2fc6015d05cbf7bcc686412a2eb100a1f2220296e3b4a644c69633a'

    archive_path = _get_archive_path(file_path, "points/bunny/bunny.npy")

    if not exists(archive_path):
        _fetch_remote(file_url, archive_path, file_checksum)
        license_path = join(split(archive_path)[0], "bunny.LICENSE")
        _fetch_remote(license_url, license_path, license_checksum)
        # Print license terms unless accept_license is set to True
        if not accept_license:
            if exists(license_path):
                with open(license_path, 'r') as f:
                    print(f.read())

    return np.load(archive_path, mmap_mode='r')
