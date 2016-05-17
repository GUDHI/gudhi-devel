ext:
	python setup.py build_ext --inplace

test:
	python test/Simplex_tree_UT.py
	python test/Rips_complex_UT.py
	python test/Alpha_complex_UT.py

example:
	python example/Simplex_tree_example.py
	python example/Rips_complex_example.py
	python example/Alpha_complex_example.py

clean:
	rm -rf build/ *.o *.so *.cpp

# For not Makefile to take into account example and test directories
.PHONY: example test
