#!/bin/csh

echo "" > out

echo -n "L1:" >> out
python rotation_analysis.py -l lens_1 -oa -j | grep "mount_position" >> out
echo -n "L2:" >> out
python rotation_analysis.py -l lens_2 -oa -j | grep "mount_position" >> out
echo -n "L3:" >> out
python rotation_analysis.py -l lens_3 -oa -j | grep "mount_position" >> out
echo -n "L4:" >> out
python rotation_analysis.py -l lens_56 -oa -j | grep "mount_position" >> out
echo -n "L56:" >> out
python rotation_analysis.py -l lens_4a -oa -j | grep "mount_position" >> out

