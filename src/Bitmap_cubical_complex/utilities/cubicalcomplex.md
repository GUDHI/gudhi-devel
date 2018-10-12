---
layout: page
title: "Cubical complex"
meta_title: "Cubical complex"
teaser: ""
permalink: /cubicalcomplex/
---
{::comment}
Leave the lines above as it is required by the web site generator 'Jekyll'
{:/comment}


## cubical_complex_persistence ##
This program computes persistent homology, by using the Bitmap_cubical_complex class, of cubical complexes provided in text files in Perseus style.
See [here](/doc/latest/fileformats.html#FileFormatsPerseus) for a description of the file format.

**Example**

```
   cubical_complex_persistence data/bitmap/CubicalTwoSphere.txt
```

* Creates a Cubical Complex from the Perseus style file `CubicalTwoSphere.txt`,
computes Persistence cohomology from it and writes the results in a persistence file `CubicalTwoSphere.txt_persistence`.

## periodic_cubical_complex_persistence ##

Same as above, but with periodic boundary conditions.

**Example**

```
   periodic_cubical_complex_persistence data/bitmap/3d_torus.txt
```

* Creates a Periodical Cubical Complex from the Perseus style file `3d_torus.txt`,
computes Persistence cohomology from it and writes the results in a persistence file `3d_torus.txt_persistence`.
