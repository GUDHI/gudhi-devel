---
layout: page
title: "Alpha complex"
meta_title: "Alpha complex"
teaser: ""
permalink: /alphacomplex/
---
{::comment}
Leave the lines above as it is required by the web site generator 'Jekyll'
{:/comment}


## alpha_complex_persistence ##

This program computes the persistent homology with coefficient field *Z/pZ* of
the dD alpha complex built from a dD point cloud.

Different versions of Alpha complex computation are available:
 * fast: right combinatorics, values can be arbitrarily bad
 * safe (default): values can have a relative error at most 1e-5
 * exact: true values rounded to double.

Default Alpha complex filtrations computation are square of the circumradius of the simplex.
If you are interested in the circumradius of the simplex as filtration values, pass the
'--squared-filtrations off' (or '-s off') option.

Alpha complex can be, or not, weighted (requires a file containing weights values).

The output diagram contains one bar per line, written with the convention:

```
   p dim birth death
```

where `dim` is the dimension of the homological feature, `birth` and `death`
are respectively the birth and death of the feature, and `p` is the
characteristic of the field *Z/pZ* used for homology coefficients (`p` must be
a prime number).

**Usage**

```
   alpha_complex_persistence [options] <input OFF file>
```

where
`<input OFF file>` is the path to the input point cloud in
[nOFF ASCII format]({{ site.officialurl }}/doc/latest/fileformats.html#FileFormatsOFF).

**Allowed options**

* `-h [ --help ]` Produce help message.
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-r [ --max-alpha-square-value ]` (default = inf) Maximal alpha square value for the Alpha complex construction.
* `-p [ --field-charac ]` (default = 11) Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value
    to see zero length intervals.
* `-w [ --weight-file ]` is the path to the file containing the weights of the points (one value per line).
    Default version is not weighted.
* `-e [ --exact ]` for the exact computation version.
* `-f [ --fast ]` for the fast computation version.
* `-s [ --squared-filtrations ]` to activate square filtration computations (default is 'on', can be 'off').

**Example**

```
   alpha_complex_persistence -p 2 -m 0.45 ../../data/points/tore3D_300.off
```

N.B.:

* Filtration values are alpha square values.
* Weights values are explained on CGAL
[dD Triangulations](https://doc.cgal.org/latest/Triangulation/index.html)
and
[Regular triangulation](https://doc.cgal.org/latest/Triangulation/index.html#TriangulationSecRT) documentation.
In this case, the filtration value of each simplex is computed as the power distance of the smallest power sphere
passing through all of its vertices. Weighted Alpha complex can have negative filtration values. This is the reason
why '-s off' or '--square-root-filtrations off' is ignored in this case (filtration values would be `NaN` in this case).

## alpha_complex_3d_persistence ##
This program computes the persistent homology with coefficient field *Z/pZ* of
the 3D Alpha complex built from a 3D point cloud.

Different versions of 3D Alpha complex computation are available:
 * fast: right combinatorics, values can be arbitrarily bad
 * safe (default): values can have a relative error at most 1e-5
 * exact: true values rounded to double.

3D Alpha complex can be, or not, weighted (requires a file containing weights values)
and/or periodic (requires a file describing the periodic domain).

The output diagram contains
one bar per line, written with the convention:

```
p dim birth death
```

where `dim` is the dimension of the homological feature, `birth` and `death`
are respectively the birth and death of the feature, and `p` is the
characteristic of the field *Z/pZ* used for homology coefficients (`p` must be
a prime number).

**Usage**

```
   alpha_complex_3d_persistence [options] <input OFF file>
```

where `<input OFF file>` is the path to the input point cloud in
[nOFF ASCII format]({{ site.officialurl }}/doc/latest/fileformats.html#FileFormatsOFF).

**Allowed options**

* `-h [ --help ]` Produce help message.
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-r [ --max-alpha-square-value ]` (default = inf) Maximal alpha square value for the Alpha complex construction.
* `-p [ --field-charac ]` (default=11) Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value
    to see zero length intervals.
* `-c [ --cuboid-file ]` is the path to the file describing the periodic domain. It must be in the format described
    [here]({{ site.officialurl }}/doc/latest/fileformats.html#FileFormatsIsoCuboid). Default version is not periodic.
* `-w [ --weight-file ]` is the path to the file containing the weights of the points (one value per line).
    Default version is not weighted.
* `-e [ --exact ]` for the exact computation version (not compatible with weight and periodic version).
* `-f [ --fast ]` for the fast computation version.

**Example**

```
alpha_complex_3d_persistence ../../data/points/tore3D_300.off -p 2 -m 0.45
```

N.B.:

* `alpha_complex_3d_persistence` only accepts OFF files in dimension 3.
* Filtration values are alpha square values.
* Weights values are explained on CGAL
[Alpha shape](https://doc.cgal.org/latest/Alpha_shapes_3/index.html#Alpha_shapes_3Definitions)
and
[Regular triangulation](https://doc.cgal.org/latest/Triangulation_3/index.html#Triangulation3secclassRegulartriangulation) documentation.
* The periodic domain is detailed on CGAL [3D Periodic Triangulations User Manual](
https://doc.cgal.org/latest/Periodic_3_triangulation_3/index.html)
