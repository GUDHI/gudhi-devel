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

def test_fetch_remote_datasets():
    # Test files download from given urls
    path_file_dw = remote.fetch("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d.csv", "spiral_2d.csv")
    names_dw = re.split(r' |/|\\', path_file_dw)
    assert 'remote_datasets' == names_dw[0]
    assert 'spiral_2d.csv' == names_dw[1]

    path_file_dw = remote.fetch("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/sphere3D_pts_on_grid.off", "sphere3D_pts_on_grid.off")
    names_dw = re.split(r' |/|\\', path_file_dw)
    assert 'remote_datasets' == names_dw[0]
    assert 'sphere3D_pts_on_grid.off' == names_dw[1]


    # Test files download with checksums provided
    path_file_dw = remote.fetch("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/spiral_2d.csv", "spiral_2d.csv", checksum_flag = True,
                                                     file_checksum = '37530355d980d957c4ec06b18c775f90a91e446107d06c6201c9b4000b077f38')
    names_dw = re.split(r' |/|\\', path_file_dw)
    assert 'remote_datasets' == names_dw[0]
    assert 'spiral_2d.csv' == names_dw[1]

    path_file_dw = remote.fetch("https://raw.githubusercontent.com/GUDHI/gudhi-data/main/points/sphere3D_pts_on_grid.off", "sphere3D_pts_on_grid.off",
                                                                      checksum_flag = True, file_checksum = '32f96d2cafb1177f0dd5e0a019b6ff5658e14a619a7815ae55ad0fc5e8bd3f88')
    names_dw = re.split(r' |/|\\', path_file_dw)
    assert 'remote_datasets' == names_dw[0]
    assert 'sphere3D_pts_on_grid.off' == names_dw[1]

    # Test spiral_2d.csv wrapping function
    path_file_dw = remote.fetch_spiral_2d()
    names_dw = re.split(r' |/|\\', path_file_dw)
    assert 'remote_datasets' == names_dw[0]
    assert 'spiral_2d.csv' == names_dw[1]
