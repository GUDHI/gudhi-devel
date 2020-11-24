# This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
# which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
# license details.
# Author(s):       Vincent Rouvreau
#
# Copyright (C) 2016 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import errno
import os
import sys

import numpy as np

from .cubical._cubical_complex import Cubical_complex_interface_, Persistent_cohomology_cubical_interface_

__author__ = "Vincent Rouvreau"
__copyright__ = "Copyright (C) 2020 Inria"
__license__ = "MIT"


# CubicalComplex python interface
class CubicalComplex:
    """The CubicalComplex is an example of a structured complex useful in
    computational mathematics (specially rigorous numerics) and image
    analysis.
    """
    thisptr = None
    pcohptr = None

    def __init__(self, dimensions=None, top_dimensional_cells=None,
                 perseus_file=''):
        """CubicalComplex constructor from dimensions and
        top_dimensional_cells or from a Perseus-style file name.

        :param dimensions: A list of number of top dimensional cells.
        :type dimensions: list of int
        :param top_dimensional_cells: A list of cells filtration values.
        :type top_dimensional_cells: list of double

        Or

        :param top_dimensional_cells: A multidimensional array of cells
            filtration values.
        :type top_dimensional_cells: anything convertible to a numpy ndarray

        Or

        :param perseus_file: A Perseus-style file name.
        :type perseus_file: string
        """
        if ((dimensions is not None) and (top_dimensional_cells is not None)
            and (perseus_file == '')):
            self.thisptr = Cubical_complex_interface_(dimensions, top_dimensional_cells)
        elif ((dimensions is None) and (top_dimensional_cells is not None)
            and (perseus_file == '')):
            top_dimensional_cells = np.array(top_dimensional_cells,
                                             copy = False,
                                             order = 'F')
            dimensions = top_dimensional_cells.shape
            top_dimensional_cells = top_dimensional_cells.ravel(order='F')
            self.thisptr = Cubical_complex_interface_(dimensions, top_dimensional_cells)
        elif ((dimensions is None) and (top_dimensional_cells is None)
            and (perseus_file != '')):
            if os.path.isfile(perseus_file):
                self.thisptr = Cubical_complex_interface_(perseus_file.encode('utf-8'))
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                        perseus_file)
        else:
            raise TypeError("CubicalComplex can be constructed from dimensions and "
                "top_dimensional_cells or from a Perseus-style file name.")

    def __dealloc__(self):
        self.thisptr = 0

    def num_simplices(self):
        """This function returns the number of all cubes in the complex.

        :returns:  int -- the number of all cubes in the complex.
        """
        return self.thisptr.num_simplices()

    def dimension(self):
        """This function returns the dimension of the complex.

        :returns:  int -- the complex dimension.
        """
        return self.thisptr.dimension()

    def compute_persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function computes the persistence of the complex, so it can be
        accessed through :func:`persistent_betti_numbers`,
        :func:`persistence_intervals_in_dimension`, etc. This function is
        equivalent to :func:`persistence` when you do not want the list
        :func:`persistence` returns.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :returns: Nothing.
        """
        self.pcohptr = Persistent_cohomology_cubical_interface_(self.thisptr, True)
        self.pcohptr.compute_persistence(homology_coeff_field, min_persistence)

    def persistence(self, homology_coeff_field=11, min_persistence=0):
        """This function computes and returns the persistence of the complex.

        :param homology_coeff_field: The homology coefficient field. Must be a
            prime number
        :type homology_coeff_field: int.
        :param min_persistence: The minimum persistence value to take into
            account (strictly greater than min_persistence). Default value is
            0.0.
            Sets min_persistence to -1.0 to see all values.
        :type min_persistence: float.
        :returns: list of pairs(dimension, pair(birth, death)) -- the
            persistence of the complex.
        """
        self.compute_persistence(homology_coeff_field, min_persistence)
        return self.pcohptr.get_persistence()

    def cofaces_of_persistence_pairs(self):
        """A persistence interval is described by a pair of cells, one that creates the
        feature and one that kills it. The filtration values of those 2 cells give coordinates
        for a point in a persistence diagram, or a bar in a barcode. Structurally, in the
        cubical complexes provided here, the filtration value of any cell is the minimum of the
        filtration values of the maximal cells that contain it. Connecting persistence diagram
        coordinates to the corresponding value in the input (i.e. the filtration values of
        the top-dimensional cells) is useful for differentiation purposes.

        This function returns a list of pairs of top-dimensional cells corresponding to
        the persistence birth and death cells of the filtration. The cells are represented by
        their indices in the input list of top-dimensional cells (and not their indices in the
        internal datastructure that includes non-maximal cells). Note that when two adjacent
        top-dimensional cells have the same filtration value, we arbitrarily return one of the two
        when calling the function on one of their common faces.

        :returns: The top-dimensional cells/cofaces of the positive and negative cells,
            together with the corresponding homological dimension, in two lists of numpy arrays of integers.
            The first list contains the regular persistence pairs, grouped by dimension.
            It contains numpy arrays of shape [number_of_persistence_points, 2].
            The indices of the arrays in the list correspond to the homological dimensions, and the
            integers of each row in each array correspond to: (index of positive top-dimensional cell,
            index of negative top-dimensional cell).
            The second list contains the essential features, grouped by dimension.
            It contains numpy arrays of shape [number_of_persistence_points, 1].
            The indices of the arrays in the list correspond to the homological dimensions, and the
            integers of each row in each array correspond to: (index of positive top-dimensional cell).
        """
        assert self.pcohptr != None, "compute_persistence() must be called before cofaces_of_persistence_pairs()"

        output = [[],[]]
        pr = np.array(self.pcohptr.cofaces_of_cubical_persistence_pairs())

        ess_ind = np.argwhere(pr[:,2] == -1)[:,0]
        ess = pr[ess_ind]
        max_h = max(ess[:,0])+1 if len(ess) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(ess[:,0] == h)[:,0]
            output[1].append(ess[hidxs][:,1])

        reg_ind = np.setdiff1d(np.array(range(len(pr))), ess_ind)
        reg = pr[reg_ind]
        max_h = max(reg[:,0])+1 if len(reg) > 0 else 0
        for h in range(max_h):
            hidxs = np.argwhere(reg[:,0] == h)[:,0]
            output[0].append(reg[hidxs][:,1:])

        return output

    def betti_numbers(self):
        """This function returns the Betti numbers of the complex.

        :returns: list of int -- The Betti numbers ([B0, B1, ..., Bn]).

        :note: betti_numbers function requires :func:`compute_persistence` function to be
            launched first.

        :note: betti_numbers function always returns [1, 0, 0, ...] as infinity
            filtration cubes are not removed from the complex.
        """
        assert self.pcohptr != None, "compute_persistence() must be called before betti_numbers()"
        return self.pcohptr.betti_numbers()

    def persistent_betti_numbers(self, from_value, to_value):
        """This function returns the persistent Betti numbers of the complex.

        :param from_value: The persistence birth limit to be added in the
            numbers (persistent birth <= from_value).
        :type from_value: float.
        :param to_value: The persistence death limit to be added in the
            numbers (persistent death > to_value).
        :type to_value: float.

        :returns: list of int -- The persistent Betti numbers ([B0, B1, ...,
            Bn]).

        :note: persistent_betti_numbers function requires :func:`compute_persistence`
            function to be launched first.
        """
        assert self.pcohptr != None, "compute_persistence() must be called before persistent_betti_numbers()"
        return self.pcohptr.persistent_betti_numbers(from_value, to_value)

    def persistence_intervals_in_dimension(self, dimension):
        """This function returns the persistence intervals of the complex in a
        specific dimension.

        :param dimension: The specific dimension.
        :type dimension: int.
        :returns: The persistence intervals.
        :rtype:  numpy array of dimension 2

        :note: intervals_in_dim function requires :func:`compute_persistence` function to be
            launched first.
        """
        assert self.pcohptr != None, "compute_persistence() must be called before persistence_intervals_in_dimension()"
        return np.array(self.pcohptr.intervals_in_dimension(dimension))
