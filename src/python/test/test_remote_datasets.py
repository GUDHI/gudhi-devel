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
import os.path
import pytest

def check_dir_file_names(path_file_dw, filename, dirname):
    assert os.path.isfile(path_file_dw)

    names_dw = re.split(r' |/|\\', path_file_dw)
    assert dirname == names_dw[0]
    assert filename == names_dw[1]

def check_fetch_output(url, filename, dirname = "remote_datasets", file_checksum = None):
    path_file_dw = remote.fetch(url, filename, dirname, file_checksum)
    check_dir_file_names(path_file_dw, filename, dirname)

def test_fetch_remote_datasets():
    # Test fetch with a wrong checksum
    with pytest.raises(OSError):
        check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d.csv", "spiral_2d.csv", file_checksum = 'XXXXXXXXXX')

    # Test files download from given urls with checksums provided
    check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d.csv", "spiral_2d.csv",
                                file_checksum = '37530355d980d957c4ec06b18c775f90a91e446107d06c6201c9b4000b077f38')

    check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/sphere3D_pts_on_grid.off", "sphere3D_pts_on_grid.off",
                                file_checksum = '32f96d2cafb1177f0dd5e0a019b6ff5658e14a619a7815ae55ad0fc5e8bd3f88')

    # Test files download from given urls without checksums
    check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d.csv", "spiral_2d.csv")

    check_fetch_output("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/sphere3D_pts_on_grid.off", "sphere3D_pts_on_grid.off")

    # Test spiral_2d.csv wrapping function
    path_file_dw = remote.fetch_spiral_2d()
    check_dir_file_names(path_file_dw, 'spiral_2d.csv', 'remote_datasets')
