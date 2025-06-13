---
layout: page
title: "OFF point set generator"
meta_title: "OFF point set generator"
teaser: ""
permalink: /pointsetgenerator/
---
{::comment}
Leave the lines above as it is required by the web site generator 'Jekyll'
{:/comment}


Generates a pointset and save it in an OFF file. Command-line is:

```
off_file_from_shape_generator on|in sphere|cube|curve|torus|klein <filename> <num_points> <dimension> [radius=1.0] [seed] ...
```

Warning: "on cube" generator is not available!

**Examples**

```
off_file_from_shape_generator on sphere onSphere.off 1000 3 15.2 42
```

* Generates an onSphere.off file with 1000 points randomized on a sphere of dimension 3, radius 15.2. The seed of the random generator has been set to 42 for reproducibility.

```
off_file_from_shape_generator in sphere inSphere.off 100 2
```

* Generates an inSphere.off file with 100 points randomized in a sphere of dimension 2 (circle) and radius 1.0 (default).

```
off_file_from_shape_generator in cube inCube.off 10000 3 5.8
```

* Generates a inCube.off file with 10000 points randomized in a cube of dimension 3 and side 5.8.
