---
layout: page
title: "Čech complex"
meta_title: "Čech complex"
teaser: ""
permalink: /cechcomplex/
---
{::comment}
Leave the lines above as it is required by the web site generator 'Jekyll'
{:/comment}


# Čech complex #

## cech_persistence ##
This program computes the persistent homology with coefficient field *Z/pZ* of
a Čech complex defined on a set of input points, using Euclidean distance.

Different versions of Delaunay-Čech complex computation are available:
 * fast: right combinatorics, values can be arbitrarily bad
 * safe (default): values can have a relative error at most 1e-5
 * exact: true values rounded to double.

The output diagram contains one bar per line, written with the convention:

`p dim birth death`

where `dim` is the dimension of the homological feature, `birth` and `death`
are respectively the birth and death of the feature, and `p` is the
characteristic of the field *Z/pZ* used for homology coefficients (`p` must be
a prime number).

**Usage**

`cech_persistence [options] <input OFF file>`

where
`<input OFF file>` is the path to the input point cloud in
[nOFF ASCII format]({{ site.officialurl }}/doc/latest/fileformats.html#FileFormatsOFF).

**Allowed options**

* `-h [ --help ]` Produce help message.
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-r [ --max-radius ]` (default = inf) Maximal radius for the Čech complex construction.
* `-d [ --cpx-dimension ]` (default = 1) Maximal dimension of the Čech complex we want to compute.
* `-p [ --field-charac ]` (default = 11) Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value
    to see zero length intervals.
* `-e [ --exact ]` for the exact computation version.
* `-f [ --fast ]` for the fast computation version.

Beware: this program may use a lot of RAM and take a lot of time if `max-radius` is set to a large value.

**Example 1 with Z/2Z coefficients**

`cech_persistence ../../data/points/tore3D_1307.off -r 0.25 -m 0.5 -d 3 -p 2`

**Example 2 with Z/3Z coefficients**

`cech_persistence ../../data/points/tore3D_1307.off -r 0.25 -m 0.5 -d 3 -p 3`


## delaunay_cech_persistence ##

This program Computes the persistent homology with coefficient field *Z/pZ*
of a Delaunay-Čech complex defined on a set of input points.

Different versions of Delaunay-Čech complex computation are available:
 * fast: right combinatorics, values can be arbitrarily bad
 * safe (default): values can have a relative error at most 1e-5
 * exact: true values rounded to double.


Default Delaunay-Čech complex filtrations computation are squared radius of the MEB.
If you are interested in the radius of the MEB as filtration values, pass the
'--squared-filtrations off' (or '-s off') option.


 The output diagram contains one bar per line, written with the
convention:

`p dim birth death`

where `dim` is the dimension of the homological feature, `birth` and `death`
are respectively the birth and death of the feature, and `p` is the
characteristic of the field *Z/pZ* used for homology coefficients (`p` must be
a prime number).

**Usage**

`delaunay_cech_persistence [options] <input OFF file>`

where
`<input OFF file>` is the path to the input point cloud in
[nOFF ASCII format]({{ site.officialurl }}/doc/latest/fileformats.html#FileFormatsOFF).

**Allowed options**

* `-h [ --help ]` Produce help message.
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-r [ --max-squared-radius ]` (default = inf)  Maximal squared length of an edge for the Delaunay-Čech complex
    construction.
* `-p [ --field-charac ]` (default = 11) Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value
    to see zero length intervals.
* `-e [ --exact ]` for the exact computation version.
* `-f [ --fast ]` for the fast computation version.
* `-s [ --squared-filtrations ]` to activate square filtration computations (default is 'on', can be 'off').

**Example 1 with Z/2Z coefficients**

`delaunay_cech_persistence ../../data/points/tore3D_1307.off -r 0.25 -m 0.5 -p 2`

**Example 2 with Z/3Z coefficients**

`delaunay_cech_persistence ../../data/points/tore3D_1307.off -r 0.25 -m 0.5 -p 3`
