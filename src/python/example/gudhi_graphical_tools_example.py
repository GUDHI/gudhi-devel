#!/usr/bin/env python

""" This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
    Author(s):       Vincent Rouvreau

    Copyright (C) 2016 Inria

    Modification(s):
      - YYYY/MM Author: Description of the modification
"""

__license__ = "MIT"


import matplotlib.pyplot as plot
from gudhi import plot_persistence_barcode, plot_persistence_diagram


print("#####################################################################")
print("Show barcode persistence example")

persistence = [
    (2, (1.0, float("inf"))),
    (1, (1.4142135623730951, float("inf"))),
    (1, (1.4142135623730951, float("inf"))),
    (0, (0.0, float("inf"))),
    (0, (0.0, 1.0)),
    (0, (0.0, 1.0)),
    (0, (0.0, 1.0)),
]
plot_persistence_barcode(persistence)
plot.show()

print("#####################################################################")
print("Show diagram persistence example")

plot_persistence_diagram(persistence)
plot.show()

print("#####################################################################")
print("Show diagram persistence example with a confidence band")

plot_persistence_diagram(persistence, band=0.2)
plot.show()

print("#####################################################################")
print("Show barcode and diagram persistence side by side example")
fig, axes = plot.subplots(nrows=1, ncols=2)
plot_persistence_barcode(persistence, axes=axes[0])
plot_persistence_diagram(persistence, axes=axes[1])
fig.suptitle("barcode versus diagram")
plot.show()
