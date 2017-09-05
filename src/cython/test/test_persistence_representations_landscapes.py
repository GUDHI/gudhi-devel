import gudhi

"""
    This file is part of the Gudhi Library. The Gudhi library
    (Geometric Understanding in Higher Dimensions) is a generic C++
    library for computational topology.

    Author(s):       Pawel Dlotko

    Copyright (C) 2017  Swansea University

    This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

    This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Pawel Dlotko"
__copyright__ = "Copyright (C) 2017 Swansea University"
__license__ = "GPL v3"

epsilon = 0.0000005;

def test_check_construction_of_landscape():    
    p = gudhi.Persistence_landscape("data/file_with_diagram",0)
    q = gudhi.Persistence_landscape
    q.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram")
    assert p == q


def test_check_construction_of_landscape_form_gudhi_style_file():
    p = gudhi.Persistence_landscape("data/persistence_file_with_four_entries_per_line", 1)
    q = gudhi.Persistence_landscape
    q.load_landscape_from_file("data/persistence_file_with_four_entries_per_line_landscape");  
    assert p == q;

def test_check_computations_of_integrals():  
    p = gudhi.Persistence_landscape("data/file_with_diagram",0)
    integral = p.compute_integral_of_landscape()
    assert fabs(integral - 2.34992) <= 0.00001


def test_check_computations_of_integrals_for_each_level_separatelly():
    diag = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram");
    p = gudhi.Persistence_landscape(diag)
    integrals_for_different_levels = [0.216432,0.204763,0.188793,0.178856,0.163142,0.155015,0.143046,0.133765,0.123531,0.117393,0.111269,0.104283,0.0941308,0.0811208,0.0679001,0.0580801,0.0489647,0.0407936,0.0342599,0.02896,0.0239881,0.0171792,0.0071511,0.00462067,0.00229033,0.000195296]
    for lv in range(0, len(integrals_for_different_levels)):
        integral = p.compute_integral_of_a_level_of_a_landscape(lv);
        assert fabs(integral - integrals_fir_different_levels[lv]) <= 0.00001

def test_check_computations_of_integrals_of_powers_of_landscape():
    diag = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram")
    p = gudhi.Persistence_landscape(diag)
    integrals_fir_different_powers = [17.1692,2.34992,0.49857,0.126405,0.0355235]
    for power in range(0,5):
        integral = p.compute_integral_of_landscape(power)
        assert fabs(integral - integrals_fir_different_powers[power]) <= 0.00005
  

def test_check_computations_of_values_on_different_points():
    diag = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram")
    p = gudhi.Persistence_landscape(diag);
    assert fabs(p.compute_value_at_a_given_point(1, 0.0)) <= 0.00001 
    assert fabs(p.compute_value_at_a_given_point(1, 0.1) - 0.0692324) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(1, 0.2) - 0.163369) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(1, 0.3) - 0.217115) <= 0.00001 
    assert fabs(p.compute_value_at_a_given_point(2, 0.0)) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(2, 0.1) - 0.0633688) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(2, 0.2) - 0.122361) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(2, 0.3) - 0.195401) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(3, 0.0)) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(3, 0.1) - 0.0455386) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(3, 0.2) - 0.0954012) <= 0.00001
    assert fabs(p.compute_value_at_a_given_point(3, 0.3) - 0.185282) <= 0.00001


def test_check_computations_of_maxima_and_norms():
    diag = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram")
    p = gudhi.Persistence_landscape(diag)
    second = gudhi.Persistence_landscape
    second.load_landscape_from_file("data/file_with_landscape_from_file_with_diagram_1")
    sum_ = gudhi.Persistence_landscape()
    sum_ = p + second;
    assert fabs(p.compute_maximum() - 0.431313) <= 0.00001
    assert fabs(p.compute_norm_of_landscape(1) - 2.34992) <= 0.00001
    assert fabs(p.compute_norm_of_landscape(2) - 0.706095) <= 0.00001
    assert fabs(p.compute_norm_of_landscape(3) - 0.501867) <= 0.00001
    assert fabs(compute_distance_of_landscapes(p, sum_, 1) - 27.9323) <= 0.00005
    assert fabs(compute_distance_of_landscapes(p, sum_, 2) - 2.35199) <= 0.00001



def test_check_default_parameters_of_distances():
    diag = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram")
    p = gudhi.Persistence_landscape(diag)
    diag1 = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram_1")
    q = gudhi.Persistence_landscape(diag1)
    dist_numeric_limit_max = p.distance(q, sys.float_info.max);
    dist_infinity = p.distance(q, sys.float_info.max);
    assert dist_numeric_limit_max == dist_infinity


def test_check_computations_of_averages():
    diag = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram")
    p = gudhi.Persistence_landscape(diag)
    diag2 = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram_1")
    q = gudhi.Persistence_landscape(diag2)
    av = gudhi.Persistence_landscape
    av.compute_average({p, q})
    template_average = Persistence_landscape
    template_averagetemplate_average.load_landscape_from_file("data/average")
    assert template_average == av


def test_check_computations_of_distances():
    diag = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram")
    p = gudhi.Persistence_landscape(diag)
    diag2 = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram_1")
    q = Persistence_landscape(diag2)
    assert fabs(p.distance(q) - 25.5824) <= 0.00005
    assert fabs(p.distance(q, 2) - 2.12636) <= 0.00001
    assert fabs(p.distance(q, sys.float_info.max) - 0.359068) <= 0.00001


def test_check_computations_of_scalar_product():
    diag = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram")
    p = gudhi.Persistence_landscape(diag)
    diag2 = read_persistence_intervals_in_one_dimension_from_file("data/file_with_diagram_1")
    q = Persistence_landscape(diag2)
    assert fabs(p.compute_scalar_product(q) - 0.754498) <= 0.00001

