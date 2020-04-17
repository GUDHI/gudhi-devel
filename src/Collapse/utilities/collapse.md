---
layout: page
title: "Collapse"
meta_title: "Edge collapse"
teaser: ""
permalink: /collapse/
---
{::comment}
Leave the lines above as it is required by the web site generator 'Jekyll'
{:/comment}


## point_cloud_edge_collapse_rips_persistence ##
This program computes the Rips graph defined on a set of input points, using Euclidean distance, and collapses edges.
This program finally computes persistent homology with coefficient field *Z/pZ* of the Rips complex built on top of these collapse edges.
The output diagram contains one bar per line, written with the convention:

`p dim birth death`

where `dim` is the dimension of the homological feature, `birth` and `death` are respectively the birth and death of the feature, and `p` is the characteristic of the field *Z/pZ* used for homology coefficients (`p` must be a prime number).

**Usage**

`point_cloud_edge_collapse_rips_persistence [options] <OFF input file>`

**Allowed options**

* `-h [ --help ]` Produce help message
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-r [ --max-edge-length ]` (default = inf) Maximal length of an edge for the Rips complex construction.
* `-d [ --cpx-dimension ]` (default = 1) Maximal dimension of the Rips complex we want to compute.
* `-p [ --field-charac ]` (default = 11)     Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value to see zero length intervals.

Beware: this program may use a lot of RAM and take a lot of time if `max-edge-length` is set to a large value.

**Example 1 with Z/2Z coefficients**

`point_cloud_edge_collapse_rips_persistence ../../data/points/tore3D_1307.off -r 0.25 -m 0.5 -d 3 -p 2`

**Example 2 with Z/3Z coefficients**

`point_cloud_edge_collapse_rips_persistence ../../data/points/tore3D_1307.off -r 0.25 -m 0.5 -d 3 -p 3`


## distance_matrix_edge_collapse_rips_persistence ##

Same as `point_cloud_edge_collapse_rips_persistence` but taking a distance matrix as input.

**Usage**

`distance_matrix_edge_collapse_rips_persistence [options] <CSV input file>`

where
`<CSV input file>` is the path to the file containing a distance matrix. Can be square or lower triangular matrix. Separator is ';'.
The code do not check if it is dealing with a distance matrix. It is the user responsibility to provide a valid input.
Please refer to data/distance_matrix/lower_triangular_distance_matrix.csv for an example of a file.

**Example**

`distance_matrix_edge_collapse_rips_persistence data/distance_matrix/full_square_distance_matrix.csv -r 15 -d 3 -p 3 -m 0`

