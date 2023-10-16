""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

from gudhi import SimplexTree
import numpy as np
import pytest

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"


def test_insertion():
    st = SimplexTree()
    assert st._is_defined() == True
    assert st._is_persistence_defined() == False
    assert st.is_empty()

    # insert test
    assert st.insert([0, 1]) == True

    assert st.dimension() == 1
    assert not st.is_empty()

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
    assert st.find([3]) == False
    assert st.find([0, 3]) == False
    assert st.find([1, 3]) == False
    assert st.find([2, 3]) == False

    # filtration test
    assert st.filtration([0, 1, 2]) == 4.0
    assert st.filtration([0, 2]) == 4.0
    assert st.filtration([1, 2]) == 4.0
    assert st.filtration([2]) == 4.0
    assert st.filtration([0, 1]) == 0.0
    assert st.filtration([0]) == 0.0
    assert st.filtration([1]) == 0.0

    # skeleton test
    assert list(st.get_skeleton(2)) == [
        ([0, 1, 2], 4.0),
        ([0, 1], 0.0),
        ([0, 2], 4.0),
        ([0], 0.0),
        ([1, 2], 4.0),
        ([1], 0.0),
        ([2], 4.0),
    ]
    assert list(st.get_skeleton(1)) == [
        ([0, 1], 0.0),
        ([0, 2], 4.0),
        ([0], 0.0),
        ([1, 2], 4.0),
        ([1], 0.0),
        ([2], 4.0),
    ]
    assert list(st.get_skeleton(0)) == [([0], 0.0), ([1], 0.0), ([2], 4.0)]

    # remove_maximal_simplex test
    assert st.get_cofaces([0, 1, 2], 1) == []
    st.remove_maximal_simplex([0, 1, 2])
    assert list(st.get_skeleton(2)) == [
        ([0, 1], 0.0),
        ([0, 2], 4.0),
        ([0], 0.0),
        ([1, 2], 4.0),
        ([1], 0.0),
        ([2], 4.0),
    ]
    assert st.find([0, 1, 2]) == False
    assert st.find([0, 1]) == True
    assert st.find([0, 2]) == True
    assert st.find([0]) == True
    assert st.find([1]) == True
    assert st.find([2]) == True

    assert st.persistence(persistence_dim_max=True) == [
        (1, (4.0, float("inf"))),
        (0, (0.0, float("inf"))),
    ]
    assert st._is_persistence_defined() == True

    assert st.betti_numbers() == [1, 1]
    assert st.persistent_betti_numbers(-0.1, 10000.0) == [0, 0]
    assert st.persistent_betti_numbers(0.0, 10000.0) == [1, 0]
    assert st.persistent_betti_numbers(3.9, 10000.0) == [1, 0]
    assert st.persistent_betti_numbers(4.0, 10000.0) == [1, 1]
    assert st.persistent_betti_numbers(9999.0, 10000.0) == [1, 1]

    st = SimplexTree()
    st.insert([1])
    st.insert([2])
    st2 = SimplexTree()
    st2.insert([1, 2])
    st2.remove_maximal_simplex([1, 2])
    assert st == st2

def test_expansion():
    st = SimplexTree()
    assert st._is_defined() == True
    assert st._is_persistence_defined() == False

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

    assert list(st.get_filtration()) == [
        ([2], 0.1),
        ([3], 0.1),
        ([2, 3], 0.1),
        ([0], 0.2),
        ([0, 2], 0.2),
        ([1], 0.3),
        ([0, 1], 0.3),
        ([1, 3], 0.4),
        ([1, 2], 0.5),
        ([5], 0.6),
        ([6], 0.6),
        ([5, 6], 0.6),
        ([4], 0.7),
        ([2, 4], 0.7),
        ([0, 3], 0.8),
        ([4, 6], 0.9),
        ([3, 6], 1.0),
    ]

    st.expansion(3)
    assert st.num_vertices() == 7
    assert st.num_simplices() == 22

    assert list(st.get_filtration()) == [
        ([2], 0.1),
        ([3], 0.1),
        ([2, 3], 0.1),
        ([0], 0.2),
        ([0, 2], 0.2),
        ([1], 0.3),
        ([0, 1], 0.3),
        ([1, 3], 0.4),
        ([1, 2], 0.5),
        ([0, 1, 2], 0.5),
        ([1, 2, 3], 0.5),
        ([5], 0.6),
        ([6], 0.6),
        ([5, 6], 0.6),
        ([4], 0.7),
        ([2, 4], 0.7),
        ([0, 3], 0.8),
        ([0, 1, 3], 0.8),
        ([0, 2, 3], 0.8),
        ([0, 1, 2, 3], 0.8),
        ([4, 6], 0.9),
        ([3, 6], 1.0),
    ]


def test_automatic_dimension():
    st = SimplexTree()
    assert st._is_defined() == True
    assert st._is_persistence_defined() == False

    # insert test
    assert st.insert([0, 1, 3], filtration=0.5) == True
    assert st.insert([0, 1, 2], filtration=1.0) == True

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
    assert st._is_defined() == True
    assert st._is_persistence_defined() == False

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
    st.insert([0, 1, 6, 7], filtration=1.0)

    assert st.make_filtration_non_decreasing() == False

    # Modify specific values to test make_filtration_non_decreasing
    st.assign_filtration([0, 1, 6, 7], 0.8)
    st.assign_filtration([0, 1, 6], 0.9)
    st.assign_filtration([0, 6], 0.6)
    st.assign_filtration([3, 4, 5], 1.2)
    st.assign_filtration([3, 4], 1.1)
    st.assign_filtration([4, 5], 1.99)

    assert st.make_filtration_non_decreasing() == True

    assert st.filtration([0, 1, 6, 7]) == 1.0
    assert st.filtration([0, 1, 6]) == 1.0
    assert st.filtration([0, 1]) == 1.0
    assert st.filtration([0]) == 1.0
    assert st.filtration([1]) == 1.0
    assert st.filtration([3, 4, 5]) == 2.0
    assert st.filtration([3, 4]) == 2.0
    assert st.filtration([4, 5]) == 2.0


def test_extend_filtration():

    # Inserted simplex:
    #      5   4
    #      o   o
    #     / \ /
    #    o   o
    #   /2\ /3
    #  o   o
    #  1   0

    st = SimplexTree()
    st.insert([0, 2])
    st.insert([1, 2])
    st.insert([0, 3])
    st.insert([2, 5])
    st.insert([3, 4])
    st.insert([3, 5])
    st.assign_filtration([0], 1.0)
    st.assign_filtration([1], 2.0)
    st.assign_filtration([2], 3.0)
    st.assign_filtration([3], 4.0)
    st.assign_filtration([4], 5.0)
    st.assign_filtration([5], 6.0)

    assert list(st.get_filtration()) == [
        ([0, 2], 0.0),
        ([1, 2], 0.0),
        ([0, 3], 0.0),
        ([3, 4], 0.0),
        ([2, 5], 0.0),
        ([3, 5], 0.0),
        ([0], 1.0),
        ([1], 2.0),
        ([2], 3.0),
        ([3], 4.0),
        ([4], 5.0),
        ([5], 6.0),
    ]

    st.extend_filtration()

    assert list(st.get_filtration()) == [
        ([6], -3.0),
        ([0], -2.0),
        ([1], -1.8),
        ([2], -1.6),
        ([0, 2], -1.6),
        ([1, 2], -1.6),
        ([3], -1.4),
        ([0, 3], -1.4),
        ([4], -1.2),
        ([3, 4], -1.2),
        ([5], -1.0),
        ([2, 5], -1.0),
        ([3, 5], -1.0),
        ([5, 6], 1.0),
        ([4, 6], 1.2),
        ([3, 6], 1.4),
        ([3, 4, 6], 1.4),
        ([3, 5, 6], 1.4),
        ([2, 6], 1.6),
        ([2, 5, 6], 1.6),
        ([1, 6], 1.8),
        ([1, 2, 6], 1.8),
        ([0, 6], 2.0),
        ([0, 2, 6], 2.0),
        ([0, 3, 6], 2.0),
    ]

    dgms = st.extended_persistence(min_persistence=-1.0)
    assert len(dgms) == 4
    # Sort by (death-birth) descending - we are only interested in those with the longest life span
    for idx in range(4):
        dgms[idx] = sorted(dgms[idx], key=lambda x: (-abs(x[1][0] - x[1][1])))

    assert dgms[0][0][1][0] == pytest.approx(2.0)
    assert dgms[0][0][1][1] == pytest.approx(3.0)
    assert dgms[1][0][1][0] == pytest.approx(5.0)
    assert dgms[1][0][1][1] == pytest.approx(4.0)
    assert dgms[2][0][1][0] == pytest.approx(1.0)
    assert dgms[2][0][1][1] == pytest.approx(6.0)
    assert dgms[3][0][1][0] == pytest.approx(6.0)
    assert dgms[3][0][1][1] == pytest.approx(1.0)


def test_simplices_iterator():
    st = SimplexTree()

    assert st.insert([0, 1, 2], filtration=4.0) == True
    assert st.insert([2, 3, 4], filtration=2.0) == True

    for simplex in st.get_simplices():
        print("simplex is: ", simplex[0])
        assert st.find(simplex[0]) == True
        print("filtration is: ", simplex[1])
        assert st.filtration(simplex[0]) == simplex[1]


def test_collapse_edges():
    st = SimplexTree()

    assert st.insert([0, 1], filtration=1.0) == True
    assert st.insert([1, 2], filtration=1.0) == True
    assert st.insert([2, 3], filtration=1.0) == True
    assert st.insert([0, 3], filtration=1.0) == True
    assert st.insert([0, 2], filtration=2.0) == True
    assert st.insert([1, 3], filtration=2.0) == True

    assert st.num_simplices() == 10

    st.collapse_edges()
    assert st.num_simplices() == 9
    assert st.find([0, 2]) == False  # [1, 3] would be fine as well
    for simplex in st.get_skeleton(0):
        assert simplex[1] == 1.0


def test_reset_filtration():
    st = SimplexTree()

    assert st.insert([0, 1, 2], 3.0) == True
    assert st.insert([0, 3], 2.0) == True
    assert st.insert([3, 4, 5], 3.0) == True
    assert st.insert([0, 1, 6, 7], 4.0) == True

    # Guaranteed by construction
    for simplex in st.get_simplices():
        assert st.filtration(simplex[0]) >= 2.0

    # dimension until 5 even if simplex tree is of dimension 3 to test the limits
    for dimension in range(5, -1, -1):
        st.reset_filtration(0.0, dimension)
        for simplex in st.get_skeleton(3):
            print(simplex)
            if len(simplex[0]) < (dimension) + 1:
                assert st.filtration(simplex[0]) >= 2.0
            else:
                assert st.filtration(simplex[0]) == 0.0


def test_boundaries_iterator():
    st = SimplexTree()

    assert st.insert([0, 1, 2, 3], filtration=1.0) == True
    assert st.insert([1, 2, 3, 4], filtration=2.0) == True

    assert list(st.get_boundaries([1, 2, 3])) == [([1, 2], 1.0), ([1, 3], 1.0), ([2, 3], 1.0)]
    assert list(st.get_boundaries([2, 3, 4])) == [([2, 3], 1.0), ([2, 4], 2.0), ([3, 4], 2.0)]
    assert list(st.get_boundaries([2])) == []

    with pytest.raises(RuntimeError):
        list(st.get_boundaries([]))

    with pytest.raises(RuntimeError):
        list(st.get_boundaries([0, 4]))  # (0, 4) does not exist

    with pytest.raises(RuntimeError):
        list(st.get_boundaries([6]))  # (6) does not exist


def test_persistence_intervals_in_dimension():
    # Here is our triangulation of a 2-torus - taken from https://dioscuri-tda.org/Paris_TDA_Tutorial_2021.html
    #   0-----3-----4-----0
    #   | \   | \   | \   | \   |
    #   |   \ |   \ |    \|   \ |
    #   1-----8-----7-----1
    #   | \   | \   | \   | \   |
    #   |   \ |   \ |   \ |   \ |
    #   2-----5-----6-----2
    #   | \   | \   | \   | \   |
    #   |   \ |   \ |   \ |   \ |
    #   0-----3-----4-----0
    st = SimplexTree()
    st.insert([0, 1, 8])
    st.insert([0, 3, 8])
    st.insert([3, 7, 8])
    st.insert([3, 4, 7])
    st.insert([1, 4, 7])
    st.insert([0, 1, 4])
    st.insert([1, 2, 5])
    st.insert([1, 5, 8])
    st.insert([5, 6, 8])
    st.insert([6, 7, 8])
    st.insert([2, 6, 7])
    st.insert([1, 2, 7])
    st.insert([0, 2, 3])
    st.insert([2, 3, 5])
    st.insert([3, 4, 5])
    st.insert([4, 5, 6])
    st.insert([0, 4, 6])
    st.insert([0, 2, 6])
    st.compute_persistence(persistence_dim_max=True)

    H0 = st.persistence_intervals_in_dimension(0)
    assert np.array_equal(H0, np.array([[0.0, float("inf")]]))
    H1 = st.persistence_intervals_in_dimension(1)
    assert np.array_equal(H1, np.array([[0.0, float("inf")], [0.0, float("inf")]]))
    H2 = st.persistence_intervals_in_dimension(2)
    assert np.array_equal(H2, np.array([[0.0, float("inf")]]))
    # Test empty case
    assert st.persistence_intervals_in_dimension(3).shape == (0, 2)


def test_equality_operator():
    st1 = SimplexTree()
    st2 = SimplexTree()

    assert st1 == st2

    st1.insert([1, 2, 3], 4.0)
    assert st1 != st2

    st2.insert([1, 2, 3], 4.0)
    assert st1 == st2


def test_simplex_tree_deep_copy():
    st = SimplexTree()
    st.insert([1, 2, 3], 0.0)
    # compute persistence only on the original
    st.compute_persistence()

    st_copy = st.copy()
    assert st_copy == st
    st_filt_list = list(st.get_filtration())

    # check persistence is not copied
    assert st._is_persistence_defined() == True
    assert st_copy._is_persistence_defined() == False

    # remove something in the copy and check the copy is included in the original
    st_copy.remove_maximal_simplex([1, 2, 3])
    a_filt_list = list(st_copy.get_filtration())
    assert len(a_filt_list) < len(st_filt_list)

    for a_splx in a_filt_list:
        assert a_splx in st_filt_list

    # test double free
    del st
    del st_copy


def test_simplex_tree_deep_copy_constructor():
    st = SimplexTree()
    st.insert([1, 2, 3], 0.0)
    # compute persistence only on the original
    st.compute_persistence()

    st_copy = SimplexTree(st)
    assert st_copy == st
    st_filt_list = list(st.get_filtration())

    # check persistence is not copied
    assert st._is_persistence_defined() == True
    assert st_copy._is_persistence_defined() == False

    # remove something in the copy and check the copy is included in the original
    st_copy.remove_maximal_simplex([1, 2, 3])
    a_filt_list = list(st_copy.get_filtration())
    assert len(a_filt_list) < len(st_filt_list)

    for a_splx in a_filt_list:
        assert a_splx in st_filt_list

    # test double free
    del st
    del st_copy


def test_simplex_tree_constructor_exception():
    with pytest.raises(TypeError):
        st = SimplexTree(other="Construction from a string shall raise an exception")


def test_create_from_array():
    a = np.array([[1, 4, 13, 6], [4, 3, 11, 5], [13, 11, 10, 12], [6, 5, 12, 2]])
    st = SimplexTree.create_from_array(a, max_filtration=5.0)
    assert list(st.get_filtration()) == [([0], 1.0), ([3], 2.0), ([1], 3.0), ([0, 1], 4.0), ([1, 3], 5.0)]
    assert SimplexTree.create_from_array([[0, 10], [10, 1]]).dimension() == 1
    assert SimplexTree.create_from_array([[0, 10], [10, 1]], max_filtration=5).dimension() == 0
    assert SimplexTree.create_from_array([[0, 10], [10, 1]], max_filtration=-5).dimension() == -1


def test_insert_edges_from_coo_matrix():
    try:
        from scipy.sparse import coo_matrix
        from scipy.spatial import cKDTree
    except ImportError:
        print("Skipping, no SciPy")
        return

    st = SimplexTree()
    st.insert([1, 2, 7], 7)
    row = np.array([2, 5, 3])
    col = np.array([1, 4, 6])
    dat = np.array([1, 2, 3])
    edges = coo_matrix((dat, (row, col)))
    st.insert_edges_from_coo_matrix(edges)
    assert list(st.get_filtration()) == [
        ([1], 1.0),
        ([2], 1.0),
        ([1, 2], 1.0),
        ([4], 2.0),
        ([5], 2.0),
        ([4, 5], 2.0),
        ([3], 3.0),
        ([6], 3.0),
        ([3, 6], 3.0),
        ([7], 7.0),
        ([1, 7], 7.0),
        ([2, 7], 7.0),
        ([1, 2, 7], 7.0),
    ]

    pts = np.random.rand(100, 2)
    tree = cKDTree(pts)
    edges = tree.sparse_distance_matrix(tree, max_distance=0.15, output_type="coo_matrix")
    st = SimplexTree()
    st.insert_edges_from_coo_matrix(edges)
    assert 100 < st.num_simplices() < 1000


def test_insert_batch():
    st = SimplexTree()
    # vertices
    st.insert_batch(np.array([[6, 1, 5]]), np.array([-5.0, 2.0, -3.0]))
    # triangles
    st.insert_batch(np.array([[2, 10], [5, 0], [6, 11]]), np.array([4.0, 0.0]))
    # edges
    st.insert_batch(np.array([[1, 5], [2, 5]]), np.array([1.0, 3.0]))

    assert list(st.get_filtration()) == [
        ([6], -5.0),
        ([5], -3.0),
        ([0], 0.0),
        ([10], 0.0),
        ([0, 10], 0.0),
        ([11], 0.0),
        ([0, 11], 0.0),
        ([10, 11], 0.0),
        ([0, 10, 11], 0.0),
        ([1], 1.0),
        ([2], 1.0),
        ([1, 2], 1.0),
        ([2, 5], 4.0),
        ([2, 6], 4.0),
        ([5, 6], 4.0),
        ([2, 5, 6], 4.0),
    ]


def test_expansion_with_blocker():
    st = SimplexTree()
    st.insert([0, 1], 0)
    st.insert([0, 2], 1)
    st.insert([0, 3], 2)
    st.insert([1, 2], 3)
    st.insert([1, 3], 4)
    st.insert([2, 3], 5)
    st.insert([2, 4], 6)
    st.insert([3, 6], 7)
    st.insert([4, 5], 8)
    st.insert([4, 6], 9)
    st.insert([5, 6], 10)
    st.insert([6], 10)

    def blocker(simplex):
        try:
            # Block all simplices that contain vertex 6
            simplex.index(6)
            print(simplex, " is blocked")
            return True
        except ValueError:
            print(simplex, " is accepted")
            st.assign_filtration(simplex, st.filtration(simplex) + 1.0)
            return False

    st.expansion_with_blocker(2, blocker)
    assert st.num_simplices() == 22
    assert st.dimension() == 2
    assert st.find([4, 5, 6]) == False
    assert st.filtration([0, 1, 2]) == 4.0
    assert st.filtration([0, 1, 3]) == 5.0
    assert st.filtration([0, 2, 3]) == 6.0
    assert st.filtration([1, 2, 3]) == 6.0

    st.expansion_with_blocker(3, blocker)
    assert st.num_simplices() == 23
    assert st.dimension() == 3
    assert st.find([4, 5, 6]) == False
    assert st.filtration([0, 1, 2]) == 4.0
    assert st.filtration([0, 1, 3]) == 5.0
    assert st.filtration([0, 2, 3]) == 6.0
    assert st.filtration([1, 2, 3]) == 6.0
    assert st.filtration([0, 1, 2, 3]) == 7.0


def test_assign_filtration_missing():
    st = SimplexTree()
    with pytest.raises(ValueError):
        st.assign_filtration([1, 2], 3)

def test_prune_above_dimension():
    st = SimplexTree()

    # insert test
    assert st.insert([0, 1, 3], filtration=0.5) == True
    assert st.insert([0, 1, 2], filtration=1.0) == True

    assert st.num_vertices() == 4
    assert st.num_simplices() == 11

    assert st.dimension() == 2
    assert st.upper_bound_dimension() == 2

    assert st.prune_above_dimension(2000) == False
    assert st.dimension() == 2
    assert st.upper_bound_dimension() == 2
    assert st.num_vertices() == 4
    assert st.num_simplices() == 11

    assert st.prune_above_dimension(1) == True
    assert st.dimension() == 1
    assert st.upper_bound_dimension() == 1
    assert st.num_vertices() == 4
    assert st.num_simplices() == 9

    assert st.prune_above_dimension(0) == True
    assert st.dimension() == 0
    assert st.upper_bound_dimension() == 0
    assert st.num_vertices() == 4
    assert st.num_simplices() == 4

    assert st.prune_above_dimension(-1) == True
    assert st.dimension() == -1
    assert st.upper_bound_dimension() == -1
    assert st.num_vertices() ==  0
    assert st.num_simplices() == 0

    assert st.prune_above_dimension(-200) == False
