"""
  This script generates faux lens decentres and tilts calculated from user-inputted 
  axis positions - its primary function is to test that the data loaded into Zemax is 
  consistent with measured data.
"""
import sys

sys.path.append("../")

import numpy as np
import json

from axis import axis

# For the following test purposes, let us define a measurement system as, looking through 
# the front of the lens, x is positive to the right, y is positive upwards and z is 
# positive propagating away from the observer.


# L1 has its axis pointing up, no decentre
#
L1_data = {}
L1 = axis((0, -50, -100), 
	  (0, 50, 100),
	  pt1_radius=None, 
	  pt2_radius=None, 
	  flip_lens=False, 
	  flip_PCS_xy_axes=False, 
	  flip_PCS_x_direction=False, 
	  flip_PCS_y_direction=False, 
	  flip_PCS_z_direction=False, 
	  z_offset=0)

L1_data['mount_position'] = 1;
L1_data['axis'] = [{'x_decentre': L1.getZemaxDecentresAndTilts()[0][0],
		   'y_decentre': L1.getZemaxDecentresAndTilts()[0][1],
		   'x_tilt': L1.getZemaxDecentresAndTilts()[1][0],
		   'y_tilt': L1.getZemaxDecentresAndTilts()[1][1],
		   'axis_type': 'OPTICAL'}]
		   
# L2 has its axis pointing down, no decentre
#
L2_data = {}
L2 = axis((0, -50, 100), 
	  (0, 50, -100),
	  pt1_radius=None, 
	  pt2_radius=None, 
	  flip_lens=False, 
	  flip_PCS_xy_axes=False, 
	  flip_PCS_x_direction=False, 
	  flip_PCS_y_direction=False, 
	  flip_PCS_z_direction=False, 
	  z_offset=0)

L2_data['mount_position'] = 1;
L2_data['axis'] = [{'x_decentre': L2.getZemaxDecentresAndTilts()[0][0],
		   'y_decentre': L2.getZemaxDecentresAndTilts()[0][1],
		   'x_tilt': L2.getZemaxDecentresAndTilts()[1][0],
		   'y_tilt': L2.getZemaxDecentresAndTilts()[1][1],
		   'axis_type': 'OPTICAL'}]

# L3 has its axis pointing to the right in x only.
#
L3_data = {}
L3 = axis((-50, 0, -100), 
	  (50, 0, 100),
	  pt1_radius=None, 
	  pt2_radius=None, 
	  flip_lens=False, 
	  flip_PCS_xy_axes=False, 
	  flip_PCS_x_direction=False, 
	  flip_PCS_y_direction=False, 
	  flip_PCS_z_direction=False, 
	  z_offset=0)

L3_data['mount_position'] = 1;
L3_data['axis'] = [{'x_decentre': L3.getZemaxDecentresAndTilts()[0][0],
		   'y_decentre': L3.getZemaxDecentresAndTilts()[0][1],
		   'x_tilt': L3.getZemaxDecentresAndTilts()[1][0],
		   'y_tilt': L3.getZemaxDecentresAndTilts()[1][1],
		   'axis_type': 'OPTICAL'}]

# L4 has its axis pointing to the left in x only.
#
L4_data = {}
L4 = axis((-50, 0, 100), 
	  (50, 0, -100),
	  pt1_radius=None, 
	  pt2_radius=None, 
	  flip_lens=False, 
	  flip_PCS_xy_axes=False, 
	  flip_PCS_x_direction=False, 
	  flip_PCS_y_direction=False, 
	  flip_PCS_z_direction=False, 
	  z_offset=0)

L4_data['mount_position'] = 1;
L4_data['axis'] = [{'x_decentre': L4.getZemaxDecentresAndTilts()[0][0],
		   'y_decentre': L4.getZemaxDecentresAndTilts()[0][1],
		   'x_tilt': L4.getZemaxDecentresAndTilts()[1][0],
		   'y_tilt': L4.getZemaxDecentresAndTilts()[1][1],
		   'axis_type': 'OPTICAL'}]

# L56 has its axis facing up and to the right, with a decentre
#
L56_data = {}
L56 = axis((-35, -50, -100), 
	   (40, 25, 100),
	   pt1_radius=None, 
	   pt2_radius=None, 
	   flip_lens=False, 
	   flip_PCS_xy_axes=False, 
	   flip_PCS_x_direction=False, 
	   flip_PCS_y_direction=False, 
	   flip_PCS_z_direction=False, 
	   z_offset=0)

L56_data['mount_position'] = 1;
L56_data['axis'] = [{'x_decentre': L56.getZemaxDecentresAndTilts()[0][0],
		    'y_decentre': L56.getZemaxDecentresAndTilts()[0][1],
		    'x_tilt': L56.getZemaxDecentresAndTilts()[1][0],
		    'y_tilt': L56.getZemaxDecentresAndTilts()[1][1],
		    'axis_type': 'OPTICAL'}]

# Output this in a format that can be loaded easily into Zemax

with open("out", 'w') as f:
  f.write("L1:[" + json.dumps(L1_data) + "]\n")
  f.write("L2:[" + json.dumps(L2_data) + "]\n")
  f.write("L3:[" + json.dumps(L3_data) + "]\n")
  f.write("L4:[" + json.dumps(L4_data) + "]\n")
  f.write("L56:[" + json.dumps(L56_data) + "]\n")






