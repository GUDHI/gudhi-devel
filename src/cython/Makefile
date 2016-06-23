ext:
	python setup.py build_ext --inplace

test:
	python test/test_simplex_tree.py
	python test/test_mini_simplex_tree.py
	python test/test_rips_complex.py
	python test/test_alpha_complex.py
	python test/test_witness_complex.py
	python -m pytest test/test_cubical_complex.py

example:
	python example/simplex_tree_example.py
	python example/mini_simplex_tree_example.py
	python example/rips_complex_example.py
	python example/alpha_complex_example.py
	python example/random_cubical_complex_example.py 10 10
	python example/cubical_complex_from_perseus_file_example.py -f ../data/bitmap/CubicalTwoSphere.txt
	python example/witness_complex_from_file_example.py data/500_random_points_on_3D_Torus.csv
	# python example/gudhi_graphical_tools_example.py
	# python example/rips_persistence_diagram.py
	# python example/alpha_complex_from_file_example.py data/500_random_points_on_3D_Torus.csv
	# python example/rips_complex_from_file_example.py data/2000_random_points_on_3D_Torus.csv

clean:
	rm -rf build/ *.o *.so *.cpp

# For not Makefile to take into account example and test directories
.PHONY: example test
