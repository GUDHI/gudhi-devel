# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

from os.path import join, exists
from os import makedirs

from urllib.request import urlretrieve
import hashlib

import numpy as np

def _checksum_sha256(file_path):
    """
    Compute the file checksum using sha256.

    Parameters
    ----------
    file_path: string
        Full path of the created file.

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

def _fetch_remote(url, filename, dirname = "remote_datasets", file_checksum = None, accept_license = False):
    """
    Fetch the wanted dataset from the given url and save it in file_path.

    Parameters
    ----------
    url : string
        The url to fetch the dataset from.
    filename : string
        The name to give to downloaded file.
    dirname : string
        The directory to save the file to. Default is "remote_datasets".
    file_checksum : string
        The file checksum using sha256 to check against the one computed on the downloaded file.
        Default is 'None'.
    accept_license : boolean
        Flag to specify if user accepts the file LICENSE and prevents from printing the corresponding license terms.
        Default is False.

    Returns
    -------
    file_path: string
        Full path of the created file.
    """

    file_path = join(dirname, filename)

    # Get the file
    urlretrieve(url, file_path)

    if file_checksum is not None:
        checksum = _checksum_sha256(file_path)
        if file_checksum != checksum:
            raise IOError("{} has a SHA256 checksum : {}, "
                        "different from expected : {}."
                        "The file may be corrupted or the given url may be wrong !".format(file_path, checksum, file_checksum))

    # Print license terms unless accept_license is set to True
    if not accept_license:
        license_file = join(dirname, "LICENSE")
        if exists(license_file) and (file_path != license_file):
            with open(license_file, 'r') as f:
                print(f.read())

    return file_path

def fetch_spiral_2d(filename = "spiral_2d.npy", dirname = "remote_datasets/spiral_2d"):
    """
    Fetch "spiral_2d.npy" remotely.

    Parameters
    ----------
    filename : string
        The name to give to downloaded file. Default is "spiral_2d.npy".
    dirname : string
        The directory to save the file to. Default is "remote_datasets/spiral_2d".

    Returns
    -------
    points: array
        Array of points stored in "spiral_2d.npy".
    """
    file_url = "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d/spiral_2d.npy"
    file_checksum = '88312ffd6df2e2cb2bde9c0e1f962d7d644c6f58dc369c7b377b298dacdc4eaf'

    archive_path = join(dirname, filename)

    if not exists(archive_path):
        # Create directory if not existing
        if not exists(dirname):
            makedirs(dirname)

        file_path_pkl = _fetch_remote(file_url, filename, dirname, file_checksum)

        return np.load(file_path_pkl, mmap_mode='r')
    else:
        return np.load(archive_path, mmap_mode='r')

def fetch_bunny(filename = "bunny.npy", dirname = "remote_datasets/bunny", accept_license = False):
    """
    Fetch "bunny.npy" remotely and its LICENSE file.

    Parameters
    ----------
    filename : string
        The name to give to downloaded file. Default is "bunny.npy".
    dirname : string
        The directory to save the file to. Default is "remote_datasets/bunny".
    accept_license : boolean
        Flag to specify if user accepts the file LICENSE and prevents from printing the corresponding license terms.
        Default is False.

    Returns
    -------
    points: array
        Array of points stored in "bunny.npy".
    """

    file_url = "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/bunny/bunny.npy"
    file_checksum = '13f7842ebb4b45370e50641ff28c88685703efa5faab14edf0bb7d113a965e1b'
    license_url = "https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/bunny/LICENSE"
    license_checksum = 'b763dbe1b2fc6015d05cbf7bcc686412a2eb100a1f2220296e3b4a644c69633a'

    archive_path = join(dirname, filename)

    if not exists(archive_path):
        # Create directory if not existing
        if not exists(dirname):
            makedirs(dirname)

        license_path = _fetch_remote(license_url, "LICENSE", dirname, license_checksum)
        file_path_pkl = _fetch_remote(file_url, filename, dirname, file_checksum, accept_license)

        return np.load(file_path_pkl, mmap_mode='r')
    else:
        return np.load(archive_path, mmap_mode='r')
