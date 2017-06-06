#!/bin/csh

echo "" > out

echo -n "1:" >> out
python rotation_analysis.py -l lens_1 -oa -j | grep "mount_position" >> out
echo -n "2:" >> out
python rotation_analysis.py -l lens_2 -oa -j | grep "mount_position" >> out
echo -n "3:" >> out
python rotation_analysis.py -l lens_3 -oa -j | grep "mount_position" >> out
echo -n "4:" >> out
python rotation_analysis.py -l lens_56 -oa -j | grep "mount_position" >> out
echo -n "56:" >> out
python rotation_analysis.py -l lens_4a -oa -j | grep "mount_position" >> out

