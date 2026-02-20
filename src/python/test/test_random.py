# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):  Vincent Rouvreau
#
# Copyright (C) 2025 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification


from gudhi.random import GudhiBitGenerator
import numpy as np
from numpy.testing import assert_almost_equal


def test_gudhi_bit_generator_default():
    assert np.random.Generator(GudhiBitGenerator()).random() != np.random.Generator(GudhiBitGenerator()).random()


def test_gudhi_bit_generator_when_seeded():
    # 42 is the seed
    first_array = np.random.Generator(GudhiBitGenerator(42)).random(5)
    second_array = np.random.Generator(GudhiBitGenerator(42)).random(5)
    print(f"{first_array} == {second_array} ?")
    assert_almost_equal(first_array, second_array, decimal=1)
