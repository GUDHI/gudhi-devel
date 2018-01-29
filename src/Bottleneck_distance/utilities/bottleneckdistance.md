---
layout: page
title: "Bottleneck distance"
meta_title: "bottleneckdistance"
subheadline: ""
teaser: ""
permalink: "/bottleneckdistance/"
---
{::comment}
These flags above are here for web site generation, please let them.
cf. https://gitlab.inria.fr/GUDHI/website
Must be in conformity with _data/navigation.yml
{:/comment}



## bottleneck_read_file_example ##

This program computes the Bottleneck distance between two persistence diagram files.

**Usage**

```
   bottleneck_read_file_example <file_1.pers> <file_2.pers> [<tolerance>]
```

where

* `<file_1.pers>` and `<file_2.pers>` must be in the format described [here](../doc/latest/fileformats.html#FileFormatsPers).
* `<tolerance>` is an error bound on the bottleneck distance (set by default to the smallest positive double value).
