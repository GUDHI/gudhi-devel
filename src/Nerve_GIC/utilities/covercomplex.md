---
layout: page
title: "Cover complex"
meta_title: "covercomplex"
subheadline: ""
teaser: ""
permalink: "/covercomplex/"
---
{::comment}
These flags above are here for web site generation, please let them.
cf. https://gitlab.inria.fr/GUDHI/website
Must be in conformity with _data/navigation.yml
{:/comment}



## Nerve ##
This program builds the Nerve of a point cloud sampled on an OFF file.
The cover C comes from the preimages of intervals covering a coordinate function,
which are then refined into their connected components using the triangulation of the .OFF file.

The program also writes a file SC.txt.
The first three lines in this file are the location of the input point cloud and the function used to compute the cover.
The fourth line contains the number of vertices nv and edges ne of the Nerve. The next nv lines represent the vertices.
Each line contains the vertex ID, the number of data points it contains, and their average color function value.
Finally, the next ne lines represent the edges, characterized by the ID of their vertices.

**Usage**

`Nerve <OFF input file> coordinate resolution gain [--v]`

where

* `coordinate` is the coordinate function to cover
* `resolution` is the number of the intervals
* `gain` is the gain for each interval
* `--v` is optional, it activates verbose mode.

**Example**

`Nerve ../../data/points/human.off 2 10 0.3`

* Builds the Nerve of a point cloud sampled on a 3D human shape (human.off).
The cover C comes from the preimages of intervals (10 intervals with gain 0.3) covering the height function (coordinate 2).

`python KeplerMapperVisuFromTxtFile.py -f ../../data/points/human.off_sc.txt`

* Constructs `human.off_sc.html` file. You can now use your favorite web browser to visualize it.

## VoronoiGIC ##

This util builds the Graph Induced Complex (GIC) of a point cloud.
It subsamples *N* points in the point cloud, which act as seeds of a geodesic VoronoÃ¯ diagram.
Each cell of the diagram is then an element of C.

The program also writes a file `*_sc.off`, that is an OFF file that can be visualized with GeomView.

**Usage**

`VoroniGIC <OFF input file> samples_number [--v]`

where

* `samples_number` is the number of samples to take from the point cloud
* `--v` is optional, it activates verbose mode.

**Example**

`VoroniGIC ../../data/points/human.off 700`

* Builds the Voronoi Graph Induced Complex with 700 subsamples from `human.off` file.
`../../data/points/human_sc.off` can be visualized with GeomView.

