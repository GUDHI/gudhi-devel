---
layout: page
title: "Reduced Vietoris-Rips degree-1 persistence"
meta_title: "Reduced Vietoris-Rips degree-1 persistence"
teaser: ""
permalink: /reducedrips/
---
{::comment}
Leave the lines above as it is required by the web site generator 'Jekyll'
{:/comment}


## reduced_rips_persistence ##
This program computes the degree-1 (i.e. *H<sub>1</sub>*) Vietoris-Rips persistent homology, with coefficients in *Z/2Z*, of either a Euclidean point cloud read from an OFF file or an arbitrary symmetric distance matrix read from a CSV file. It uses the *Reduced Vietoris-Rips filtration*, which never builds the full Vietoris-Rips complex and so scales to much larger inputs than a full-complex computation would.

The output diagram contains one bar per line, written with the convention:

`p dim birth death`

where `p` is the characteristic of the coefficient field (always `2` here, as coefficients are fixed to *Z/2Z*), `dim` is the degree of the homological feature (always `1` here), and `birth` and `death` are respectively its birth and death values.

**Usage**

`reduced_rips_persistence [options] <OFF input file>`

or

`reduced_rips_persistence [options] -d <distance matrix CSV>`

**Allowed options**

* `-h [ --help ]` Produce help message.
* `-d [ --distance-matrix ]` Read a lower-triangular distance matrix (';'-separated CSV: line *i* holds the distances to points *0..i-1*). The matrix may be any symmetric, non-negative dissimilarity.
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-k [ --num-neighbors ]` (default = 0) Initial neighbor budget per point used to seed the computation. `0` means use `sqrt(n)`, which is a sensible default; the budget grows automatically when a point needs more neighbors, so this only tunes the starting allocation.
* `-s [ --search ]` (default = auto) Spatial search strategy for the point-cloud input, one of `auto | kd | brute`. `kd` uses a kd-tree for neighbor queries (best in low dimension); `brute` uses O(n^2) brute-force search (can win in higher dimension, where kd-trees degrade); `auto` picks per the ambient dimension of the cloud. Ignored when a distance matrix is given (matrix queries always scan rows).
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime (`death - birth`) of a bar to be recorded.

**Output**

In addition to the diagram, the program prints to standard error the wall-clock time and the size of the computation (number of 1-simplices, 2-simplices, and persistent pairs).

**Limitations**

* Only homological degree 1 is computed.

**Example with a point cloud**

`reduced_rips_persistence ../../data/points/tore3D_300.off -k 0`

**Example writing the diagram to a file with brute-force search**

`reduced_rips_persistence ../../data/points/tore3D_300.off -s brute -o tore3D_300_h1.pers`

**Example with a distance matrix**

`reduced_rips_persistence -d distances.csv -o diagram_h1.pers`
