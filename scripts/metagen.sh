#!/bin/bash
sep="_"
for geom in "sphere" "klein" "torus"
do
   for number in 10 100 1000
   do
      for dim in {3..5}
      do
         echo "./off_file_from_shape_generator on $geom $geom$sep$number$sep$dim.off $number $dim"
         ./off_file_from_shape_generator on $geom $geom$sep$number$sep$dim.off $number $dim
      done
   done
done

#./off_file_from_shape_generator in|on sphere|cube off_file_name points_number[integer > 0] dimension[integer > 1] radius[double > 0.0 | default = 1.0]
