#!/usr/bin/env python

import gudhi

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

print("#####################################################################")
print("Show barcode persistence example")

persistence = [(2, (1.0, float('inf'))), (1, (1.4142135623730951, float('inf'))),
               (1, (1.4142135623730951, float('inf'))), (0, (0.0, float('inf'))),
               (0, (0.0, 1.0)), (0, (0.0, 1.0)), (0, (0.0, 1.0))]
gudhi.plot_persistence_barcode(persistence)

print("#####################################################################")
print("Show diagram persistence example")

pplot = gudhi.plot_persistence_diagram(persistence)
pplot.show()

print("#####################################################################")
print("Show diagram persistence example with a confidence band")

pplot = gudhi.plot_persistence_diagram(persistence, band=0.2)
pplot.show()
