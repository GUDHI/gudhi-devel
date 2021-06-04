# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Hind Montassif
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import hashlib

from os.path import join, exists
from os import makedirs

from urllib.request import urlretrieve


def _checksum_sha256(file_path):
    """
    Compute the file checksum using sha256

    Parameters
    ----------
    file_path: string
        Full path of the created file.

    Returns
    -------
        The hex digest of file_path
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

def fetch(url, filename, dirname = "remote_datasets", file_checksum = None):
    """
    Fetch the wanted dataset from the given url and save it in file_path

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

    Returns
    -------
    file_path: string
        Full path of the created file.
    """

    file_path = join(dirname, filename)

    # Check that an existing file does not already exist at file_path
    if not exists(file_path):
        # Create directory if not existing
        if not exists(dirname):
            makedirs(dirname)

        urlretrieve(url, file_path)

        if file_checksum is not None:
            checksum = _checksum_sha256(file_path)
            if file_checksum != checksum:
                raise IOError("{} has a SHA256 checksum : {}, "
                            "different from expected : {}."
                            "The file may be corrupted or the given url may be wrong !".format(file_path, checksum,
                                                        file_checksum))

    return file_path

def fetch_spiral_2d(filename = "spiral_2d.csv", dirname = "remote_datasets"):
    """
    Fetch spiral_2d.csv remotely

    Parameters
    ----------
    filename : string
        The name to give to downloaded file. Default is "spiral_2d.csv"
    dirname : string
        The directory to save the file to. Default is "remote_datasets".

    Returns
    -------
    file_path: string
        Full path of the created file.
    """
    return fetch("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d.csv", filename, dirname,
                 '37530355d980d957c4ec06b18c775f90a91e446107d06c6201c9b4000b077f38')
