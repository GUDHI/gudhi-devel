/*! \page fileformats File formats
 
 \tableofcontents

 \section FileFormatsPers Persistence Diagram

 Such a file, whose extension is usually `.pers`, contains a list of persistence intervals.<br>
 Lines starting with `#` are ignored (comments).<br>
 Other lines might contain 2, 3 or 4 values (the number of values on each line must be the same for all lines):
 \code{.unparsed}
   [[field] dimension] birth death
 \endcode

 Here is a simple sample file:
 \code{.unparsed}
   # Beautiful persistence diagram
   2 2.7 3.7
   2 9.6 14.
   3 34.2 34.974
   4 3. inf
 \endcode

 Other sample files can be found in the `data/persistence_diagram` folder.
 
*/