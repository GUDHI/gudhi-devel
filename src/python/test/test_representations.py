import os
import sys
import matplotlib.pyplot as plt
# Disable graphics for testing purposes
plt.show = lambda:None
here = os.path.dirname(os.path.realpath(__file__))
sys.path.append(here + "/../example")
import diagram_vectorizations_distances_kernels
# pytest is unhappy if there are 0 tests
def test_nothing():
    return None
