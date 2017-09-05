import gudhi

"""This file is part of the Gudhi Library. The Gudhi library
   (Geometric Understanding in Higher Dimensions) is a generic C++
   library for computational topology.

   Author(s):       Pawel Dlotko

   Copyright (C) 2017 Swansea University

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
 
def test_check_min_max_function():	
	p = gudhi.PersistenceIntervals(None,None,"data/file_with_diagram")
	min_max_ = p.get_x_range()
	assert fabs(min_max_[0] - 0.0290362) <= epsilon
	assert fabs(min_max_[1] - 0.994537) <= epsilon

def test_check_length_of_dominant_intervals():
	p = gudhi.PersistenceIntervals(None,None,"data/file_with_diagram")
	dominant_ten_intervals_length = p.length_of_dominant_intervals(10)
	dominant_intervals_length_ = [0.862625,0.800893,0.762061,0.756501,0.729367,0.718177,0.708395,0.702844,0.700468,0.622177]	
	assert dominant_ten_intervals_length == dominant_intervals_length_ 


def test_check_dominant_intervals(): 
	p = gudhi.PersistenceIntervals(None,None,"data/file_with_diagram");
	ten_dominant_intervals = p.dominant_intervals(10);
	templ = [ (0.114718, 0.977343) , (0.133638, 0.93453) , (0.104599, 0.866659) , (0.149798, 0.906299), (0.247352, 0.976719) , (0.192675, 0.910852) , (0.191836, 0.900231) , (0.284998, 0.987842) , (0.294069, 0.994537), (0.267421, 0.889597)]	
	assert fabs(ten_dominant_intervals - templ) <= epsilon
	
	
def test_check_histogram_of_lengths():
	p = gudhi.PersistenceIntervals(None,None,"data/file_with_diagram")
	histogram = p.histogram_of_lengths(10);
	template_histogram = [10,5,3,4,4,3,6,1,7,1,1]	
	assert fabs(histogram - template_histogram) <= epsilon


def test_check_cumulative_histograms_of_lengths():
	p = gudhi.PersistenceIntervals(None,None,"data/file_with_diagram")
	cumulative_histogram = p.cumulative_histogram_of_lengths(10)
	template_cumulative_histogram  = [10,15,18,22,26,29,35,36,43,44,45]
	assert fabs(cumulative_histogram - template_cumulative_histogram) <= epsilon


def test_check_characteristic_function_of_diagram():
	p = gudhi.PersistenceIntervals(None,None,"data/file_with_diagram")
	min_max_ = p.get_x_range();
	char_funct_diag = p.characteristic_function_of_diagram(min_max_[0], min_max_[1]);
	template_char_funct_diag = [0.370665,0.84058,1.24649,1.3664,1.34032,1.31904,1.14076,0.991259,0.800714,0.0676303]	
	assert fabs(char_funct_diag - template_char_funct_diag) <= 0.0001


def test_check_cumulative_characteristic_function_of_diagram():
	p = gudhi.PersistenceIntervals(None,None,"data/file_with_diagram")
	min_max_ = p.get_x_range()
	cumul_char_funct_diag = p.cumulative_characteristic_function_of_diagram(min_max_.first, min_max_.second,None);
	template_char_funct_diag_cumul = [0.370665,1.21125,2.45774,3.82414,5.16446,6.4835,7.62426,8.61552,9.41623,9.48386]	
	assert fabs(cumul_char_funct_diag - template_char_funct_diag_cumul) <= 0.0001


def test_check_compute_persistent_betti_numbers():
	p = gudhi.PersistenceIntervals(None,None,"data/file_with_diagram")
	pbns = [(0.0290362, 1),(0.0307676, 2),(0.0366312, 3),(0.0544614, 4),(0.0920033, 5),(0.104599, 6),(0.114718, 7),(0.117379, 8),(0.123493, 9),(0.133638, 10)(0.137798, 9),(0.149798, 10),(0.155421, 11),(0.158443, 12)(0.176956, 13),(0.183234, 12),(0.191069, 13),(0.191333, 14),(0.191836, 15),(0.192675, 16),(0.208564, 17),(0.218425, 18),(0.219902, 17),(0.23233, 16),(0.234558, 17),(0.237166, 16),(0.247352, 17),(0.267421, 18),(0.268093, 19),(0.278734, 18),(0.284722, 19),(0.284998, 20),(0.294069, 21),(0.306293, 22),(0.322361, 21),(0.323152, 22),(0.371021, 23),(0.372395, 24),(0.387744, 25),(0.435537, 26),(0.462911, 25),(0.483569, 26),(0.489209, 25),(0.517115, 24),(0.522197, 23),(0.532665, 22),(0.545262, 23),(0.587227, 22),(0.593036, 23),(0.602647, 24),(0.605044, 25),(0.621962, 24),(0.629449, 23),(0.636719, 22),(0.64957, 21),(0.650781, 22),(0.654951, 23),(0.683489, 24),(0.687172, 23),(0.69703, 22),(0.701174, 21),(0.717623, 22),(0.722023, 21),(0.722298, 20),(0.725347, 19),(0.73071, 18),(0.758355, 17),(0.770913, 18),(0.790833, 17),(0.821211, 16),(0.849305, 17),(0.853669, 16),(0.866659, 15),(0.872896, 16),(0.889597, 15),(0.900231, 14),(0.903847, 13),(0.906299, 12),(0.910852, 11),(0.93453, 10),(0.944757, 9),(0.947812, 8),(0.959154, 7),(0.975654, 6),(0.976719, 5),(0.977343, 4),(0.980129, 3),(0.987842, 2),(0.990127, 1),(0.994537, 0)]
	pbns_new = p.compute_persistent_betti_numbers();  
	assert fabs(pbns - pbns_new) <= epsilon
    

