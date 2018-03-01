#!/bin/csh
#
# This shell scripts invokes rotation_analysis.py for each lens configuration, 
# producing JSON output that can be easily installed into the camera optimisation 
# package.
#

set CONFIG_FILE = "../config.json"
set OUT_FILE = "out"

echo "" > $OUT_FILE

echo "Processing Lens 1..."
echo -n "L1:" >> $OUT_FILE
python ../../rotation_analysis.py -c $CONFIG_FILE -l lens_1_redo -oa -j | grep "mount_position" >> $OUT_FILE

echo "Processing Lens 2..."
echo -n "L2:" >> $OUT_FILE
python ../../rotation_analysis.py -c $CONFIG_FILE -l lens_2_redo -oa -j | grep "mount_position" >> $OUT_FILE

echo "Processing Lens 3..."
echo -n "L3:" >> $OUT_FILE
python ../../rotation_analysis.py -c $CONFIG_FILE -l lens_3_redo -oa -j | grep "mount_position" >> $OUT_FILE

echo "Processing Lens 4..."
echo -n "L4:" >> $OUT_FILE
python ../../rotation_analysis.py -c $CONFIG_FILE -l lens_4b_redo -oa -j | grep "mount_position" >> $OUT_FILE

echo "Processing Lens 56..."
echo -n "L56:" >> $OUT_FILE
python ../../rotation_analysis.py -c $CONFIG_FILE -l lens_56_redo -oa -j | grep "mount_position" >> $OUT_FILE

