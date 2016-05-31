ext:
	python setup.py build_ext --inplace

test:
	python test/test_simplex_tree.py
	python test/test_mini_simplex_tree.py
	python test/test_rips_complex.py
	python test/test_alpha_complex.py

example:
	python example/simplex_tree_example.py
	python example/mini_simplex_tree_example.py
	python example/rips_complex_example.py
	python example/alpha_complex_example.py
	python example/cubical_complex_example.py

clean:
	rm -rf build/ *.o *.so *.cpp

# For not Makefile to take into account example and test directories
.PHONY: example test
