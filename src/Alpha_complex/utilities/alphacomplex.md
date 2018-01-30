

# Alpha complex #


## alpha_complex_persistence ##

This program computes the persistent homology with coefficient field Z/pZ of the dD alpha complex built from a dD point cloud.
The output diagram contains one bar per line, written with the convention:

```
   p dim birth death
```

where `dim` is the dimension of the homological feature, `birth` and `death` are respectively the birth and death of the feature,
and `p` is the characteristic of the field *Z/pZ* used for homology coefficients (`p` must be a prime number).

**Usage**

```
   alpha_complex_persistence [options] <input OFF file>
```

where
`<input OFF file>` is the path to the input point cloud in [nOFF ASCII format](http://www.geomview.org/docs/html/OFF.html).

**Allowed options**

* `-h [ --help ]` Produce help message
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-r [ --max-alpha-square-value ]` (default = inf) Maximal alpha square value for the Alpha complex construction.
* `-p [ --field-charac ]` (default = 11)     Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value to see zero length intervals.

**Example**

```
   alpha_complex_persistence -r 32 -p 2 -m 0.45 ../../data/points/tore3D_300.off
```

N.B.:

* Filtration values are alpha square values.


## alpha_complex_3d_persistence ##
This program computes the persistent homology with coefficient field Z/pZ of the 3D alpha complex built from a 3D point cloud. The output diagram contains one bar per line, written with the convention:

```
p dim birth death
```

where `dim` is the dimension of the homological feature, `birth` and `death` are respectively the birth and death of the feature, and `p` is the characteristic of the field *Z/pZ* used for homology coefficients (`p` must be a prime number).

**Usage**

```
   alpha_complex_3d_persistence [options] <input OFF file>
```

where `<input OFF file>` is the path to the input point cloud in [nOFF ASCII format](http://www.geomview.org/docs/html/OFF.html).

**Allowed options**

* `-h [ --help ]` Produce help message
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-p [ --field-charac ]` (default=11) Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value to see zero length intervals.

**Example**

```
alpha_complex_3d_persistence ../../data/points/tore3D_300.off -p 2 -m 0.45
```

N.B.:

* `alpha_complex_3d_persistence` only accepts OFF files in dimension 3.
* Filtration values are alpha square values.


## exact_alpha_complex_3d_persistence ##

Same as `alpha_complex_3d_persistence`, but using exact computation.
It is slower, but it is necessary when points are on a grid for instance.



## weighted_alpha_complex_3d_persistence ##

Same as `alpha_complex_3d_persistence`, but using weighted points.

**Usage**

```
   weighted_alpha_complex_3d_persistence [options] <input OFF file> <weights input file>
```

where

* `<input OFF file>` is the path to the input point cloud in [nOFF ASCII format](http://www.geomview.org/docs/html/OFF.html).
* `<input weights file>` is the path to the file containing the weights of the points (one value per line).

**Allowed options**

* `-h [ --help ]` Produce help message
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-p [ --field-charac ]` (default=11) Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value to see zero length intervals.

**Example**

```
   weighted_alpha_complex_3d_persistence ../../data/points/tore3D_300.off ../../data/points/tore3D_300.weights -p 2 -m 0.45
```


N.B.:

* Weights values are explained on CGAL [Alpha shape](https://doc.cgal.org/latest/Alpha_shapes_3/index.html#title0)
and [Regular triangulation](https://doc.cgal.org/latest/Triangulation_3/index.html#Triangulation3secclassRegulartriangulation) documentation.
* Filtration values are alpha square values.


## periodic_alpha_complex_3d_persistence ##
Same as `alpha_complex_3d_persistence`, but using periodic alpha shape 3d.
Refer to the [CGAL's 3D Periodic Triangulations User Manual](https://doc.cgal.org/latest/Periodic_3_triangulation_3/index.html) for more details.

**Usage**

```
   periodic_alpha_complex_3d_persistence [options] <input OFF file> <cuboid file>
```

where

* `<input OFF file>` is the path to the input point cloud in [nOFF ASCII format](http://www.geomview.org/docs/html/OFF.html).
* `<cuboid file>` is the path to the file describing the periodic domain. It must be in the format described
[here](/doc/latest/fileformats.html#FileFormatsIsoCuboid).

**Allowed options**

* `-h [ --help ]` Produce help message
* `-o [ --output-file ]` Name of file in which the persistence diagram is written. Default print in standard output.
* `-p [ --field-charac ]` (default=11) Characteristic p of the coefficient field Z/pZ for computing homology.
* `-m [ --min-persistence ]` (default = 0) Minimal lifetime of homology feature to be recorded. Enter a negative value to see zero length intervals


**Example**

```
periodic_alpha_complex_3d_persistence ../../data/points/grid_10_10_10_in_0_1.off ../../data/points/iso_cuboid_3_in_0_1.txt -p 3 -m 1.0
```

N.B.:

* Cuboid file must be in the format described [here](/doc/latest/fileformats.html#FileFormatsIsoCuboid).
* Filtration values are alpha square values.
