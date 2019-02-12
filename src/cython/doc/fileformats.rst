:orphan:

.. To get rid of WARNING: document isn't included in any toctree

File formats
############

Persistence Diagram
*******************

Such a file, whose extension is usually ``.pers``, contains a list of
persistence intervals.

Lines starting with ``#`` are ignored (comments).

Other lines might contain 2, 3 or 4 values (the number of values on each line
must be the same for all lines)::

    [[field] dimension] birth death

Here is a simple sample file::

    # Persistence diagram example
    2 2.7 3.7
    2 9.6 14.
    # Some comments
    3 34.2 34.974
    4 3. inf

Other sample files can be found in the `data/persistence_diagram` folder.

Such files can be generated with
:meth:`gudhi.SimplexTree.write_persistence_diagram`, read with
:meth:`gudhi.read_persistence_intervals_grouped_by_dimension`, or
:meth:`gudhi.read_persistence_intervals_in_dimension` and displayed with
:meth:`gudhi.plot_persistence_barcode` or
:meth:`gudhi.plot_persistence_diagram`.

Iso-cuboid
**********

Such a file describes an iso-oriented cuboid with diagonal opposite vertices
(min_x, min_y, min_z,...) and (max_x, max_y, max_z, ...). The format is::

    min_x min_y [min_z ...]
    max_x max_y [max_z ...]

Here is a simple sample file in the 3D case::

    -1. -1. -1.
    1. 1. 1.


.. _Perseus file format:

Perseus
*******

This file format is a format inspired from the
`Perseus software <http://www.sas.upenn.edu/~vnanda/perseus/>`_ by Vidit Nanda.
The first line contains a number d begin the dimension of the bitmap (2 in the
example below). Next d lines are the numbers of top dimensional cubes in each
dimensions (3 and 3 in the example below). Next, in lexicographical order, the
filtration of top dimensional cubes is given (1 4 6 8 20 4 7 6 5 in the example
below).

.. figure::
    ../../doc/Bitmap_cubical_complex/exampleBitmap.png
    :alt: Example of a input data.
    :figclass: align-center

    Example of a input data.

The input file for the following complex is:

.. literalinclude:: ../../data/bitmap/cubicalcomplexdoc.txt

.. centered:: ../../data/bitmap/cubicalcomplexdoc.txt

To indicate periodic boundary conditions in a given direction, then number of
top dimensional cells in this direction have to be multiplied by -1. For
instance:

.. literalinclude:: ../../data/bitmap/periodiccubicalcomplexdoc.txt

.. centered:: ../../data/bitmap/periodiccubicalcomplexdoc.txt


Indicate that we have imposed periodic boundary conditions in the direction x,
but not in the direction y.

Other sample files can be found in the `data/bitmap` folder.

.. note::
    Unlike in Perseus format the filtration on the maximal cubes can be any
    double precision number. Consequently one cannot mark the cubes that are
    not present with ``-1``'s. To do that please set their filtration value to
    :math:`+\infty` (aka. ``inf`` in the file).