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

This program computes the persistent homology with coefficient field Z/pZ of
the dD alpha complex built from a dD point cloud.
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
[nOFF ASCII format](http://www.geomview.org/docs/html/OFF.html).

**Allowed options**

* `-h [ --help ]` Produce help message
* `-o [ --output-file ]` Name of file in which the persistence diagram is
written. Default print in standard output.
* `-r [ --max-alpha-square-value ]` (default = inf) Maximal alpha square value
for the Alpha complex construction.
* `-p [ --field-charac ]` (default = 11)     Characteristic p of the
coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature
to be recorded. Enter a negative value to see zero length intervals.

**Example**

```
   alpha_complex_persistence -r 32 -p 2 -m 0.45 ../../data/points/tore3D_300.off
```

N.B.:

* Filtration values are alpha square values.


## alpha_complex_3d_persistence ##
This program computes the persistent homology with coefficient field Z/pZ of
the 3D alpha complex built from a 3D point cloud.
One can use exact computation. It is slower, but it is necessary when points
are on a grid for instance.
Alpha complex 3d can be weighted and/or periodic (refer to the
[CGAL's 3D Periodic Triangulations User Manual](
https://doc.cgal.org/latest/Periodic_3_triangulation_3/index.html)
for more details).

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
[nOFF ASCII format](http://www.geomview.org/docs/html/OFF.html).

**Allowed options**

* `-h [ --help ]` Produce help message
* `-o [ --output-file ]` Name of file in which the persistence diagram is
written. Default print in standard output.
* `-r [ --max-alpha-square-value ]` (default = inf) Maximal alpha square value
for the Alpha complex construction.
* `-p [ --field-charac ]` (default=11) Characteristic p of the coefficient
field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature
to be recorded. Enter a negative value to see zero length intervals.
* `-c [ --cuboid-file ]` is the path to the file describing the periodic domain.
It must be in the format described
[here]({{ site.officialurl }}/doc/latest/fileformats.html#FileFormatsIsoCuboid).
* `-w [ --weight-file ]` is the path to the file containing the weights of the
points (one value per line).
* `-e [ --exact ]` for the exact computation version (not compatible with
weight and periodic version).

**Example**

```
alpha_complex_3d_persistence ../../data/points/tore3D_300.off -p 2 -m 0.45
```

N.B.:

* `alpha_complex_3d_persistence` only accepts OFF files in dimension 3.
* Filtration values are alpha square values.
* Weights values are explained on CGAL
[Alpha shape](https://doc.cgal.org/latest/Alpha_shapes_3/index.html#title0)
and
[Regular triangulation](https://doc.cgal.org/latest/Triangulation_3/index.html#Triangulation3secclassRegulartriangulation) documentation.
