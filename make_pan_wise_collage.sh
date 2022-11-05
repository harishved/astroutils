#!/bin/bash
# MAKE A COLLAGE OF PANSTARRS AND WISE CUT OUTS AT A GIVEN RA AND DEC POSITION
cat 1.txt | while read line; do 
   echo $line
   name=`python3 get_name.py $line`;
   echo $name;
   mkdir -p $name
   python3 get_wise.py $line
   rm -rf unwise*fits $name/
   mv ${name}*png $name/
   python3 get_panstarrs.py $line
   rm -rf cutout*fits $name/
   mv ${name}*png $name/
   montage ${name}/${name}_g.png ${name}/${name}_r.png ${name}/${name}_i.png ${name}/${name}_z.png ${name}/${name}_y.png ${name}/${name}_W1.png ${name}/${name}_W2.png -mode Concatenate -tile x1 ${name}_montage.png
done
