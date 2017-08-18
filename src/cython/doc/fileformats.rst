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

Other sample files can be found in the data/persistence_diagram folder.

Such files can be generated with
:meth:`gudhi.SimplexTree.write_persistence_diagram`, read with
:meth:`gudhi.read_persistence_intervals_grouped_by_dimension`, or
:meth:`gudhi.read_persistence_intervals_in_dimension` and displayed with
:meth:`gudhi.plot_persistence_barcode` or
:meth:`gudhi.plot_persistence_diagram`.
