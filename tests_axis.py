import numpy as np

from axis import axis

# consistency of input and angular values

ax1 = axis((0, -50, -50), 
	    (0, 50, 50),
	    pt1_radius=None, 
	    pt2_radius=None, 
	    flip_lens=False, 
	    flip_PCS_xy_axes=False, 
	    flip_PCS_x_direction=False, 
	    flip_PCS_y_direction=False, 
	    flip_PCS_z_direction=False, 
	    z_offset=0)

ax1_decentres_and_tilts = ax1.getZemaxDecentresAndTilts()

ax2 = axis((-50, 0, -50), 
	    (50, 0, 50),
	    pt1_radius=None, 
	    pt2_radius=None, 
	    flip_lens=False, 
	    flip_PCS_xy_axes=True, 
	    flip_PCS_x_direction=False, 
	    flip_PCS_y_direction=False, 
	    flip_PCS_z_direction=False, 
	    z_offset=0)

ax2_decentres_and_tilts = ax2.getZemaxDecentresAndTilts()

ax3 = axis((-50, 0, 50), 
	    (50, 0, -50),
	    pt1_radius=None, 
	    pt2_radius=None, 
	    flip_lens=False, 
	    flip_PCS_xy_axes=True, 
	    flip_PCS_x_direction=False, 
	    flip_PCS_y_direction=False, 
	    flip_PCS_z_direction=True, 
	    z_offset=0)

ax3_decentres_and_tilts = ax3.getZemaxDecentresAndTilts()

assert ax1_decentres_and_tilts[0] == ax2_decentres_and_tilts[0] == ax3_decentres_and_tilts[0]
assert tuple(ax1_decentres_and_tilts[1]) == tuple(ax2_decentres_and_tilts[1]) == tuple(ax3_decentres_and_tilts[1])
assert np.allclose(tuple(ax1_decentres_and_tilts[1]), (45, 0, 0))

# decentre values

ax4 = axis((-50, -50, -50), 
	    (50, 50, 100),
	    pt1_radius=None, 
	    pt2_radius=None, 
	    flip_lens=False, 
	    flip_PCS_xy_axes=True, 
	    flip_PCS_x_direction=False, 
	    flip_PCS_y_direction=False, 
	    flip_PCS_z_direction=True, 
	    z_offset=0)

ax4_decentres_and_tilts = ax4.getZemaxDecentresAndTilts()

assert np.allclose(tuple(ax4_decentres_and_tilts[0]), (-16.66666, -16.66666))

