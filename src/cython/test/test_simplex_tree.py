from gudhi import SimplexTree

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_insertion():
    st = SimplexTree()
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    # insert test
    assert st.insert([0, 1]) == True

    assert st.dimension() == 1

    assert st.insert([0, 1, 2], filtration=4.0) == True

    assert st.dimension() == 2

    assert st.num_simplices() == 7
    assert st.num_vertices() == 3

    # find test
    assert st.find([0, 1, 2]) == True
    assert st.find([0, 1]) == True
    assert st.find([0, 2]) == True
    assert st.find([0]) == True
    assert st.find([1]) == True
    assert st.find([2]) == True
    assert st.find([3])    == False
    assert st.find([0, 3]) == False
    assert st.find([1, 3]) == False
    assert st.find([2, 3]) == False

    # filtration test
    st.initialize_filtration()
    assert st.filtration([0, 1, 2]) == 4.0
    assert st.filtration([0, 2]) == 4.0
    assert st.filtration([1, 2]) == 4.0
    assert st.filtration([2]) == 4.0
    assert st.filtration([0, 1]) == 0.0
    assert st.filtration([0]) == 0.0
    assert st.filtration([1]) == 0.0

    # skeleton test
    assert st.get_skeleton(2) == \
        [([0, 1, 2], 4.0), ([0, 1], 0.0), ([0, 2], 4.0),
        ([0], 0.0), ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)]
    assert st.get_skeleton(1) == \
        [([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0),
        ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)]
    assert st.get_skeleton(0) == \
        [([0], 0.0), ([1], 0.0), ([2], 4.0)]

    # remove_maximal_simplex test
    assert st.get_cofaces([0, 1, 2], 1) == []
    st.remove_maximal_simplex([0, 1, 2])
    assert st.get_skeleton(2) == \
        [([0, 1], 0.0), ([0, 2], 4.0), ([0], 0.0),
        ([1, 2], 4.0), ([1], 0.0), ([2], 4.0)]
    assert st.find([0, 1, 2]) == False
    assert st.find([0, 1]) == True
    assert st.find([0, 2]) == True
    assert st.find([0]) == True
    assert st.find([1]) == True
    assert st.find([2]) == True

    st.initialize_filtration()
    assert st.persistence(persistence_dim_max = True) == [(1, (4.0, float('inf'))), (0, (0.0, float('inf')))]
    assert st.__is_persistence_defined() == True

    assert st.betti_numbers() == [1, 1]
    assert st.persistent_betti_numbers(-0.1, 10000.0) == [0, 0]
    assert st.persistent_betti_numbers(0.0, 10000.0) == [1, 0]
    assert st.persistent_betti_numbers(3.9, 10000.0) == [1, 0]
    assert st.persistent_betti_numbers(4.0, 10000.0) == [1, 1]
    assert st.persistent_betti_numbers(9999.0, 10000.0) == [1, 1]

def test_expansion():
    st = SimplexTree()
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    # insert test
    assert st.insert([3, 2], 0.1) == True
    assert st.insert([2, 0], 0.2) == True
    assert st.insert([1, 0], 0.3) == True
    assert st.insert([3, 1], 0.4) == True
    assert st.insert([2, 1], 0.5) == True
    assert st.insert([6, 5], 0.6) == True
    assert st.insert([4, 2], 0.7) == True
    assert st.insert([3, 0], 0.8) == True
    assert st.insert([6, 4], 0.9) == True
    assert st.insert([6, 3], 1.0) == True

    assert st.num_vertices() == 7
    assert st.num_simplices() == 17
    assert st.get_filtration() == [([2], 0.1), ([3], 0.1), ([2, 3], 0.1),
    ([0], 0.2), ([0, 2], 0.2), ([1], 0.3), ([0, 1], 0.3), ([1, 3], 0.4),
    ([1, 2], 0.5), ([5], 0.6), ([6], 0.6), ([5, 6], 0.6), ([4], 0.7),
    ([2, 4], 0.7), ([0, 3], 0.8), ([4, 6], 0.9), ([3, 6], 1.0)]

    st.expansion(3)
    assert st.num_vertices() == 7
    assert st.num_simplices() == 22
    st.initialize_filtration()

    assert st.get_filtration() == [([2], 0.1), ([3], 0.1), ([2, 3], 0.1),
    ([0], 0.2), ([0, 2], 0.2), ([1], 0.3), ([0, 1], 0.3), ([1, 3], 0.4),
    ([1, 2], 0.5), ([0, 1, 2], 0.5), ([1, 2, 3], 0.5), ([5], 0.6), ([6], 0.6),
    ([5, 6], 0.6), ([4], 0.7), ([2, 4], 0.7), ([0, 3], 0.8), ([0, 1, 3], 0.8),
    ([0, 2, 3], 0.8), ([0, 1, 2, 3], 0.8), ([4, 6], 0.9), ([3, 6], 1.0)]

def test_automatic_dimension():
    st = SimplexTree()
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    # insert test
    assert st.insert([0,1,3], filtration=0.5) == True
    assert st.insert([0,1,2], filtration=1.) == True

    assert st.num_vertices() == 4
    assert st.num_simplices() == 11

    assert st.dimension() == 2
    assert st.upper_bound_dimension() == 2

    assert st.prune_above_filtration(0.6) == True
    assert st.dimension() == 2
    assert st.upper_bound_dimension() == 2

    st.assign_filtration([0, 1, 3], 0.7)
    assert st.filtration([0, 1, 3]) == 0.7

    st.remove_maximal_simplex([0, 1, 3])
    assert st.upper_bound_dimension() == 2
    assert st.dimension() == 1
    assert st.upper_bound_dimension() == 1

def test_make_filtration_non_decreasing():
    st = SimplexTree()
    assert st.__is_defined() == True
    assert st.__is_persistence_defined() == False

    # Inserted simplex:
    #    1
    #    o
    #   /X\
    #  o---o---o---o
    #  2   0   3\X/4
    #            o
    #            5
    assert st.insert([2, 1, 0], filtration=2.0) == True
    assert st.insert([3, 0], filtration=2.0) == True
    assert st.insert([3, 4, 5], filtration=2.0) == True

    assert st.make_filtration_non_decreasing() == False

    # Because of non decreasing property of simplex tree, { 0 } , { 1 } and
    # { 0, 1 } are going to be set from value 2.0 to 1.0
    st.insert([0, 1, 6, 7], filtration=1.0);

    assert st.make_filtration_non_decreasing() == False

    # Modify specific values to test make_filtration_non_decreasing
    st.assign_filtration([0,1,6,7], 0.8);
    st.assign_filtration([0,1,6], 0.9);
    st.assign_filtration([0,6], 0.6);
    st.assign_filtration([3,4,5], 1.2);
    st.assign_filtration([3,4], 1.1);
    st.assign_filtration([4,5], 1.99);

    assert st.make_filtration_non_decreasing() == True

    assert st.filtration([0,1,6,7]) == 1.
    assert st.filtration([0,1,6]) == 1.
    assert st.filtration([0,1]) == 1.
    assert st.filtration([0]) == 1.
    assert st.filtration([1]) == 1.
    assert st.filtration([3,4,5]) == 2.
    assert st.filtration([3,4]) == 2.
    assert st.filtration([4,5]) == 2.
