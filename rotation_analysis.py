#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import argparse
import datetime
from decimal import Decimal
import json
import pprint

import numpy as np
from tabulate import tabulate

from database import CMM_access
from errors import measurementError as m_err
from plotter import sag as p_sag
from rotation import sag as rot_sag, hysteresis as rot_hys
from axis import axis

def go(args, cfg):
  # Find information regarding this configuration in the config file.
  #
  configuration = [c for c in cfg['CONFIGURATIONS'] if c['id'] == args.l]
  try:
    assert len(configuration) == 1
    configuration = configuration[0]
  except:
    print "Either no configuration, or multiple configurations have been found with the same id."
    exit(0)
  rotation_data = [d for d in cfg['DATA'] if d['id'] == configuration['rotation_data']]
  error_data = [d for d in cfg['DATA'] if d['id'] == configuration['error_data']]

  # Retrieve the corresponding data for rotation/error calculations.
  #
  try:
    assert len(rotation_data) == 1
    assert len(error_data) == 1
    rotation_data = rotation_data[0]
    error_data = error_data[0]
  except:
    print "Either the datasets specified in the configuration doesn't exist, "\
      "or multiple datasets have been found with the same id."
    exit(0)
    
  try:
    db = CMM_access(cfg['PARAMS']['db_path'])
  except TypeError:
    print "TypeError excepted. Most likely reason is that the database does not exist at the "\
      "location specifed in the configuration file."
  
  #
  # SINGLE POSITION ERRORS
  # ----------------------
  #
  # Go through each entry in the error data (multiple measurements of single rotation positions) and work 
  # out the PCS transformed XY position of the optical and mechanical axes at z=0. We then use this as a 
  # measure of the error in both measuring and fitting.
  #
  OA_err_axes = []
  MA_err_axes = []
  err_mount_angles = []
  err_datetimes = []
  this_elMsRecNr_range = error_data['elMsRecNr_range']
  for elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
    for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
      element = db.getElementFromelRecNr(elRecNr)
      if element['elId'] == configuration['lens_front_elId']:		# front lens
	LENS_FRONT_CENTRE_XYZ = db.transElActPos1IntoPCS(element['elRecNr'], configuration['error_data_csId'])
	LENS_FRONT_RADIUS = db.getElementFromelRecNr(element['elRecNr'])['elActDim1']
	err_datetimes.append(db.getElementFromelRecNr(element['elRecNr'])['elTime'])
      elif element['elId'] == configuration['lens_rear_elId']:		# rear lens
	LENS_REAR_CENTRE_XYZ = db.transElActPos1IntoPCS(element['elRecNr'], configuration['error_data_csId'])
	LENS_REAR_RADIUS = db.getElementFromelRecNr(element['elRecNr'])['elActDim1']
      elif element['elId'] == configuration['mount_front_elId']:	# front mount
	MNT_FRONT_XYZ = db.transElActPos1IntoPCS(element['elRecNr'], configuration['error_data_csId'])
      elif element['elId'] == configuration['mount_rear_elId']:		# rear mount
	MNT_REAR_XYZ = db.transElActPos1IntoPCS(element['elRecNr'], configuration['error_data_csId'])
	
    if args.oa:
      # optical axis
      try:
        OA_err_axes.append(axis(LENS_FRONT_CENTRE_XYZ, LENS_REAR_CENTRE_XYZ, 
				LENS_FRONT_RADIUS, LENS_REAR_RADIUS,
				flip_lens=configuration['flip_lens'],
				flip_PCS_z_direction=configuration['flip_PCS_z_direction'],
				mount_ring_thickness=configuration['mount_ring_thickness']))
        
      except IndexError:
	print "Optical axis error data is empty."
	exit(0)
      except UnboundLocalError:
	print "Optical axis error data is empty."
	exit(0)

    if args.ma:
      # mechanical axis
      try:
        MA_err_axes.append(axis(MNT_FRONT_XYZ, MNT_REAR_XYZ, 
				None, None,
				flip_lens=configuration['flip_lens'],
				flip_PCS_z_direction=configuration['flip_PCS_z_direction'],
				mount_ring_thickness=configuration['mount_ring_thickness']))
      except IndexError:
	print "Mechanical axis error data is empty."
	exit(0)
      except UnboundLocalError:
	print "Mechanical axis error data is empty."
	exit(0)
	
    # Get the angle by matching the elMsRecNr with the corresponding entry in the COORDINATE_SYSTEMS 
    # section of the configuration file.
    #
    csys = [c for c in cfg['COORDINATE_SYSTEMS'] if c['csMsRecNr'] == elMsRecNr]
    try:
      assert len(csys) == 1
    except:
      print "Either couldn't find a coordinate system with a csMsRecNr of " + str(elMsRecNr) + ", "\
        "or multiple coordinate systems have been found."
    mount_angle_to_position_index = {
      0: 1,
      60: 2,
      120: 3,
      180: 4,
      240: 5,
      300: 6,
      360: 1
      }
    err_mount_angles.append(mount_angle_to_position_index[csys[0]['angle']])	

  # Calculation of optical axis metrics:
  #
  #  sag - xy position of axis at the z-centre of the mount ring for each orientation.
  #  angular deviation from [0, 0, 1] - angle between OA and z-axis.
  #  hysteresis - difference between starting and ending xy position evaluated at the z-centre 
  #               of the mount ring for each orientation.
  #  12 o'clock angle - angular rotation of the vector in the XY plane drawn from the fit centre 
  #                     to the xy position (evaluated at z-centre of mount ring), relative to the 
  #                     reference vector [0, 1] (clockwise from 12 o'clock).
  #  tait-bryan angles - individual rotations required to align the OA vector with the 
  #                      coordinate axis (order is ZYX)
  #  thickness - lens thickness
  #  radii - lens radii
  #  
  # and error (err*) calculations.
  #
  if args.oa:
    OA_err_xy_zisMidMountRing = []
    OA_err_angles_from_mount_normal = []
    OA_err_xy_angles_from_12_o_clock = [] 
    OA_err_xy_euler_angles = []
    OA_err_lens_thicknesses = []
    OA_err_lens_front_radii = []
    OA_err_lens_rear_radii = []
    for ax in OA_err_axes:
      OA_err_xy_zisMidMountRing.append((ax.getXY()))
      OA_err_angles_from_mount_normal.append(ax.getComponentAngles(inDeg=True))
      OA_err_xy_euler_angles.append(ax.getEulerAngles(align_axis=[0, 0, 1]))
      OA_err_lens_thicknesses.append(ax.getLensCentreThickness())
      OA_err_lens_front_radii.append(ax.pt1_radius)
      OA_err_lens_rear_radii.append(ax.pt2_radius)
  
    err = m_err(OA_err_axes) 
    OA_err_x_y_error, OA_r_err_euclidean_error = err.calculate_position_error_at_z()
    OA_err_angle_from_mount_normal_error = err.calculate_angle_error_at_z()
    OA_err_lens_thickness_error = err.calculate_lens_thickness_error()
    OA_err_lens_front_radii_error = err.calculate_lens_radius1_error()
    OA_err_lens_rear_radii_error = err.calculate_lens_radius2_error()
  
  # Calculation of mechanical axis metrics:
  #
  #  sag - xy position of axis at the z-centre of the mount ring for each orientation.
  #  angular deviation from [0, 0, 1] - angle between OA and z-axis.
  #  hysteresis - difference between starting and ending xy position evaluated at the z-centre 
  #               of the mount ring for each orientation.
  #  12 o'clock angle - angular rotation of the vector in the XY plane drawn from the fit centre 
  #                     to the xy position (evaluated at z-centre of mount ring), relative to the 
  #                     reference vector [0, 1] (clockwise from 12 o'clock).
  #  
  # and error (err*) calculations.
  #
  MA_err_xy_zisMidMountRing = []
  MA_err_angles_from_mount_normal = []
  MA_err_xy_angles_from_12_o_clock = []
  MA_err_xy_euler_angles = []
  if args.ma:
    for ax in MA_err_axes:
      MA_err_xy_zisMidMountRing.append((ax.getXY()))
      MA_err_angles_from_mount_normal.append(ax.getComponentAngles(inDeg=True))
      MA_err_xy_euler_angles.append(ax.getEulerAngles(align_axis=[0, 0, 1]))

    err = m_err(MA_err_axes) 
    MA_err_x_y_error, MA_r_err_euclidean_error = err.calculate_position_error_at_z() 
    MA_err_angle_from_mount_normal_error = err.calculate_angle_error_at_z()
  
  #
  # ANALYSIS
  # --------
  #
  # Go through each entry in the rotation data and work out the PCS transformed XY position of the 
  # optical and mechanical axes.
  #
  optical_axes = []
  mechanical_axes = []
  mount_angles = []
  datetimes = []
  this_elMsRecNr_range = rotation_data['elMsRecNr_range']
  for elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
    for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
      element = db.getElementFromelRecNr(elRecNr)
      if element['elId'] == configuration['lens_front_elId']:		# front lens
	LENS_FRONT_CENTRE_XYZ = db.transElActPos1IntoPCS(element['elRecNr'], configuration['rotation_data_csId'])
	LENS_FRONT_RADIUS = db.getElementFromelRecNr(element['elRecNr'])['elActDim1']
	datetimes.append(db.getElementFromelRecNr(element['elRecNr'])['elTime'])
      elif element['elId'] == configuration['lens_rear_elId']:		# rear lens
	LENS_REAR_CENTRE_XYZ = db.transElActPos1IntoPCS(element['elRecNr'], configuration['rotation_data_csId'])
	LENS_REAR_RADIUS = db.getElementFromelRecNr(element['elRecNr'])['elActDim1']
      elif element['elId'] == configuration['mount_front_elId']:	# front mount
	MNT_FRONT_XYZ = db.transElActPos1IntoPCS(element['elRecNr'], configuration['rotation_data_csId'])
      elif element['elId'] == configuration['mount_rear_elId']:		# rear mount
	MNT_REAR_XYZ = db.transElActPos1IntoPCS(element['elRecNr'], configuration['rotation_data_csId'])
     
    # Optical axis	
    if args.oa:
      try:
        OA = axis(LENS_FRONT_CENTRE_XYZ, LENS_REAR_CENTRE_XYZ, 
		  LENS_FRONT_RADIUS, LENS_REAR_RADIUS,
		  flip_lens=configuration['flip_lens'],
		  flip_PCS_z_direction=configuration['flip_PCS_z_direction'],
		  mount_ring_thickness=configuration['mount_ring_thickness'])	
      except IndexError:
	print "Optical axis error data is empty."
	exit(0)
      except UnboundLocalError:
	print "Optical axis error data is empty."
	exit(0)
      optical_axes.append(OA)
    
    # Mechanical axis
    if args.ma:
      try:
        MA = axis(MNT_FRONT_XYZ, MNT_REAR_XYZ, 
		  None, None,
		  flip_lens=configuration['flip_lens'],
		  flip_PCS_z_direction=configuration['flip_PCS_z_direction'],
		  mount_ring_thickness=configuration['mount_ring_thickness'])
      except IndexError:
	print "Mechanical axis error data is empty."
	exit(0)
      except UnboundLocalError:
	print "Mechanical axis error data is empty."
	exit(0)
      mechanical_axes.append(MA)
      
    # Get the angle by matching the elMsRecNr with the corresponding entry in the COORDINATE_SYSTEMS 
    # section of the configuration file.
    #
    csys = [c for c in cfg['COORDINATE_SYSTEMS'] if c['csMsRecNr'] == elMsRecNr]
    try:
      assert len(csys) == 1
    except:
      print "Either couldn't find a coordinate system with a csMsRecNr of " + str(elMsRecNr) + ", "\
        "or multiple coordinate systems have been found."
    mount_angle_to_position_index = {
      0: 1,
      60: 2,
      120: 3,
      180: 4,
      240: 5,
      300: 6,
      360: 1
      }
    mount_angles.append(mount_angle_to_position_index[csys[0]['angle']])
   
  # Calculation of optical axis metrics:
  #
  #  sag - xy position of axis at the z-centre of the mount ring for each orientation.
  #  angular deviation from [0, 0, 1] - angle between OA and z-axis.
  #  hysteresis - difference between starting and ending xy position evaluated at the z-centre 
  #               of the mount ring for each orientation.
  #  12 o'clock angle - angular rotation of the vector in the XY plane drawn from the fit centre 
  #                     to the xy position (evaluated at z-centre of mount ring), relative to the 
  #                     reference vector [0, 1] (clockwise from 12 o'clock).
  #  tait-bryan angles - individual rotations required to align the OA vector with the 
  #                      coordinate axis (order is ZYX)
  #  thickness - lens thickness
  #  radii - lens radii
  #
  OA_xy_zisMidMountRing = []
  OA_angles_from_mount_normal = []
  OA_xy_angles_from_12_o_clock = [] 
  OA_xy_euler_angles = []
  OA_lens_thicknesses = []
  OA_lens_front_radii = []
  OA_lens_rear_radii = []
  if args.oa:
    sag = rot_sag(optical_axes)
    OA_r_sag = sag.calculate()
    
    hys = rot_hys(optical_axes, mount_angles)
    OA_r_hys = hys.calculate(configuration['hys_idx_1'], configuration['hys_idx_2'])  
    
    for ax in optical_axes:
      OA_xy_zisMidMountRing.append((ax.getXY()))
      OA_angles_from_mount_normal.append(ax.getComponentAngles(inDeg=True))
      OA_xy_angles_from_12_o_clock.append(ax.getProjectedAngleInXYPlane(ref_axis=[0,1], centre=OA_r_sag[0:2]))
      OA_xy_euler_angles.append(ax.getEulerAngles(align_axis=[0, 0, 1]))
      OA_lens_thicknesses.append(ax.getLensCentreThickness())
      OA_lens_front_radii.append(ax.pt1_radius)
      OA_lens_rear_radii.append(ax.pt2_radius)
  
  # Calculation of mechanical axis metrics:
  #
  #  sag - xy position of axis at the z-centre of the mount ring for each orientation.
  #  angular deviation from [0, 0, 1] - angle between OA and z-axis.
  #  hysteresis - difference between starting and ending xy position evaluated at the z-centre 
  #               of the mount ring for each orientation.
  #  12 o'clock angle - angular rotation of the vector in the XY plane drawn from the fit centre 
  #                     to the xy position (evaluated at z-centre of mount ring), relative to the 
  #                     reference vector [0, 1] (clockwise from 12 o'clock).
  MA_xy_zisMidMountRing = []
  MA_angles_from_mount_normal = []
  MA_xy_angles_from_12_o_clock = []
  MA_xy_euler_angles = []
  if args.ma:
    sag = rot_sag(mechanical_axes)
    MA_r_sag = sag.calculate()
    
    hys = rot_hys(mechanical_axes, mount_angles)
    MA_r_hys = hys.calculate(configuration['hys_idx_1'], configuration['hys_idx_2'])  
    
    for ax in mechanical_axes:
      MA_xy_zisMidMountRing.append((ax.getXY()))
      MA_angles_from_mount_normal.append(ax.getComponentAngles(inDeg=True))
      MA_xy_angles_from_12_o_clock.append(ax.getProjectedAngleInXYPlane(ref_axis=[0,1], centre=MA_r_sag[0:2]))
      MA_xy_euler_angles.append(ax.getEulerAngles(align_axis=[0, 0, 1]))

  # Print any information, if requested.
  if args.p or args.r or args.j:
    headers1 = ['AXIS',
		'TYPE',
		'DATETIME',
	        'ORIENTATION',
	        'X_C (micron)', 
	        'Y_C (micron)', 
	        'Z AXIS ANGLE (\')',
	        'XY ANGLE (deg)',
	        'THICKNESS (mm)',
	        'FRONT RAD. (mm)',
	        'REAR RAD. (mm)',
	        'X TILT (\')',
	        'Y TILT (\')'
	       ]
    
    # per mount orientation
    #
    data1 = []
    if args.oa:
      for (datetime,
	   mount_angle, 
	   xy, 
	   mount_normal_angle,
	   thickness, 
	   front_radius, 
	   rear_radius, 
	   euler_angle) in zip(err_datetimes,
			    err_mount_angles, 
			    OA_err_xy_zisMidMountRing, 
			    OA_err_angles_from_mount_normal,
			    OA_err_lens_thicknesses, 
			    OA_err_lens_front_radii,
			    OA_err_lens_rear_radii,
			    OA_err_xy_euler_angles):
	data1.append(['OPTICAL',
	            'S',
	            datetime.strftime("%Y-%m-%d %H:%M:%S"),
		    round(mount_angle),
		    str(round(xy[0]*10**3, 1)) + " +/- " + str(round(OA_err_x_y_error[0]*10**3, 1)),
		    str(round(xy[1]*10**3, 1)) + " +/- " + str(round(OA_err_x_y_error[1]*10**3, 1)),
		    str(round(mount_normal_angle[2]*60, 2)) + " +/- " + str(round(OA_err_angle_from_mount_normal_error*60, 2)),
		    "N/A",
		    str(round(thickness, 3)) + " +/- " + str(round(OA_err_lens_thickness_error*10**3, 1)),
		    str(round(front_radius, 3)) + " +/- " + str(round(OA_err_lens_front_radii_error*10**3, 1)),
		    str(round(rear_radius, 3)) + " +/- " + str(round(OA_err_lens_rear_radii_error*10**3, 1)),
		    str(round(np.degrees(euler_angle[0])*60, 1)),
	            str(round(np.degrees(euler_angle[1])*60, 1))
		    ])	      
      for (datetime,
	   mount_angle, 
	   xy, 
	   mount_normal_angle, 
	   xy_angle, 
	   thickness, 
	   front_radius, 
	   rear_radius, 
	   euler_angle) in zip(datetimes,
			    mount_angles, 
			    OA_xy_zisMidMountRing, 
			    OA_angles_from_mount_normal, 
			    OA_xy_angles_from_12_o_clock, 
			    OA_lens_thicknesses, 
			    OA_lens_front_radii,
			    OA_lens_rear_radii,
			    OA_xy_euler_angles):
	data1.append(['OPTICAL',
	            'R',
	            datetime.strftime("%Y-%m-%d %H:%M:%S"),
		    round(mount_angle),
		    str(round(xy[0]*10**3, 1)) + " +/- " + str(round(OA_err_x_y_error[0]*10**3, 1)),
		    str(round(xy[1]*10**3, 1)) + " +/- " + str(round(OA_err_x_y_error[1]*10**3, 1)),
		    str(round(mount_normal_angle[2]*60, 2)) + " +/- " + str(round(OA_err_angle_from_mount_normal_error*60, 2)),
		    round(xy_angle),
		    str(round(thickness, 3)) + " +/- " + str(round(OA_err_lens_thickness_error*10**3, 1)),
		    str(round(front_radius, 3)) + " +/- " + str(round(OA_err_lens_front_radii_error*10**3, 1)),
		    str(round(rear_radius, 3)) + " +/- " + str(round(OA_err_lens_rear_radii_error*10**3, 1)),
		    str(round(np.degrees(euler_angle[0])*60, 1)),
	            str(round(np.degrees(euler_angle[1])*60, 1))
		    ])
    if args.ma:
      for (datetime,
	   mount_angle, 
	   xy, 
	   mount_normal_angle,
	   euler_angle) in zip(err_datetimes,
			    err_mount_angles, 
			    MA_err_xy_zisMidMountRing, 
			    MA_err_angles_from_mount_normal,
			    MA_err_xy_euler_angles):
	data1.append(['MECHANICAL',
	            'S',
	            datetime.strftime("%Y-%m-%d %H:%M:%S"),
		    round(mount_angle),
		    str(round(xy[0]*10**3, 1)) + " +/- " + str(round(MA_err_x_y_error[0]*10**3, 1)),
		    str(round(xy[1]*10**3, 1)) + " +/- " + str(round(MA_err_x_y_error[1]*10**3, 1)),
		    str(round(mount_normal_angle[2]*60, 2)) + " +/- " + str(round(MA_err_angle_from_mount_normal_error*60, 2)),
		    "N/A",
		    "N/A",
		    "N/A",
		    "N/A",
		    str(round(np.degrees(euler_angle[0])*60, 1)),
	            str(round(np.degrees(euler_angle[1])*60, 1))
		    ]) 
      for (datetime,
	   mount_angle, 
	   xy, 
	   mount_normal_angle,
	   xy_angle,
	   euler_angle) in zip(datetimes,
			    mount_angles, 
			    MA_xy_zisMidMountRing, 
			    MA_angles_from_mount_normal, 
			    MA_xy_angles_from_12_o_clock,
			    MA_xy_euler_angles):
	data1.append(['MECHANICAL',
	            'R',
	            datetime.strftime("%Y-%m-%d %H:%M:%S"),
		    round(mount_angle),
		    str(round(xy[0]*10**3, 1)) + " +/- " + str(round(MA_err_x_y_error[0]*10**3, 1)),
		    str(round(xy[1]*10**3, 1)) + " +/- " + str(round(MA_err_x_y_error[1]*10**3, 1)),
		    str(round(mount_normal_angle[2]*60, 2)) + " +/- " + str(round(MA_err_angle_from_mount_normal_error*60, 2)),
		    round(xy_angle),
		    "N/A",
		    "N/A",
		    "N/A",
		    str(round(np.degrees(euler_angle[0])*60, 1)),
	            str(round(np.degrees(euler_angle[1])*60, 1))
		    ]) 	

    # axes
    #
    headers2 = ['AXIS LABEL',
	       'FIT XY CENTRE (micron)', 
	       'FIT RADIUS (micron)', 
	       'FIT RESIDUAL RMS (micron)',   
	       'HYSTERESIS (micron)', 
	       ]
    data2 = []
    if args.oa:
      data2.append(['OPTICAL',
	      tuple((round(OA_r_sag[0]*10**3, 1), round(OA_r_sag[1]*10**3, 1))), 
	      str(round(OA_r_sag[2]*10**3, 1)), 
	      str(round(np.mean(((OA_r_sag[3]*10**3)**2))**0.5,1)),
	      str(round(OA_r_hys*10**3, 1)), 
	      ])
    if args.ma:
      data2.append(['MECHANICAL',
	      tuple((round(MA_r_sag[0]*10**3, 1), round(MA_r_sag[1]*10**3, 1))), 
	      str(round(MA_r_sag[2]*10**3, 1)), 
	      str(round(np.mean(((MA_r_sag[3]*10**3)**2))**0.5,1)),
	      str(round(MA_r_hys*10**3, 1)), 
	      ])
       
  # Now we can plot, if requested. We construct datasets first in case we want to plot optical and 
  # mechanical results on the same axes.
  #
  if args.p2d:
    datasets = []
    if args.oa:
      datasets.append({'heading': 'Optical Axis',
		      'data': {
			'x': [xy[0] for xy in OA_xy_zisMidMountRing],
			'y': [xy[1] for xy in OA_xy_zisMidMountRing],
			'x_err': np.amax(OA_err_x_y_error),
			'y_err': np.amax(OA_err_x_y_error),
			'fit_xc': OA_r_sag[0],
			'fit_yc': OA_r_sag[1],
			'fit_r': OA_r_sag[2],
			'mount_angles': np.array(np.radians(mount_angles)),
			'residuals': OA_r_sag[3],
			'angles_from_mount_normal': OA_angles_from_mount_normal,
			'xy_angles_from_12_o_clock': np.array(np.radians(OA_xy_angles_from_12_o_clock))
			
			}
		      })
    if args.ma:
      datasets.append({'heading': 'Mechanical Axis',
		      'data': {
			'x': [xy[0] for xy in MA_xy_zisMidMountRing],
			'y': [xy[1] for xy in MA_xy_zisMidMountRing],
			'x_err': np.amax(MA_err_x_y_error),
			'y_err': np.amax(MA_err_x_y_error),
			'fit_xc': MA_r_sag[0],
			'fit_yc': MA_r_sag[1],
			'fit_r': MA_r_sag[2],
			'mount_angles': np.array(np.radians(mount_angles)),
			'residuals': MA_r_sag[3],
			'angles_from_mount_normal': MA_angles_from_mount_normal,
			'xy_angles_from_12_o_clock': np.array(np.radians(MA_xy_angles_from_12_o_clock))
			}
		      })
    p = p_sag(datasets)
    p.plot()
  
  if args.p:
    print
    print "*************"
    print "PRETTY FORMAT"
    print "*************"
    print
    print "[" + configuration['id'] + "]"
    print 
    print tabulate(data1, headers1, numalign='left', stralign='left')     
    print 
    print tabulate(data2, headers2, numalign='left', stralign='left')     
    print
    
  if args.p2d:
    p.draw(hard=False)
    
  if args.r:
    print
    print "*************"
    print "REPORT FORMAT"
    print "*************"
    print
    print '\t'.join(headers1)
    for r in data1:
      print '\t'.join([str(v) for v in r])
    print
    print '\t'.join(headers2)
    for r in data2:
      print '\t'.join([str(v) for v in r])
    print
      
  if args.j:
    print
    print "***********"
    print "JSON FORMAT"
    print "***********"
    print 
    data = []
    for key, value in mount_angle_to_position_index.iteritems():
      insert = True
      for entry in data:
	if value == entry['mount_position']:
	  insert = False
      if insert:
	data.append({"mount_position": value, "axis" : []})

    for entry in data1[:-1:]:
      if entry[1] == 'R':
	for json_entry in data:
	  if json_entry['mount_position'] == entry[3]:
	    axis_type = entry[0]
	    x_decentre = "{:.2E}".format(Decimal(entry[4].split()[0].strip())/(10**3))
	    y_decentre = "{:.2E}".format(Decimal(entry[5].split()[0].strip())/(10**3))
	    x_tilt = "{:.2E}".format(Decimal(entry[11].split()[0].strip())/60)
	    y_tilt = "{:.2E}".format(Decimal(entry[12].split()[0].strip())/60)
	    insert = True;
	    for existing_axis in json_entry['axis']:
	      if axis_type == existing_axis['axis_type']:
		insert = False
	    if insert: 
	      json_entry['axis'].append({"axis_type": axis_type, 
				  "x_decentre": x_decentre,
				  "y_decentre": y_decentre,
				  "x_tilt": x_tilt,
				  "y_tilt": y_tilt})
    print json.dumps(data)

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-l", help="configuration id to process", default="lens_1")
  parser.add_argument("-p", help="pretty print output data?", action='store_true')
  parser.add_argument("-r", help="print output data in report format?", action='store_true')
  parser.add_argument("-j", help="print output data in .json format?", action="store_true")
  parser.add_argument("-p2d", help="plot 2d?", action='store_true')
  parser.add_argument("-oa", help="consider optical axis?", action='store_true')
  parser.add_argument("-ma", help="consider mechanical axis?", action='store_true')
  parser.add_argument("-d", help="debug?", action='store_true')
  parser.add_argument("-c", help="configuration file", action='store', default='config.json')
  args = parser.parse_args()
 
  with open(args.c) as fp:
    cfg = json.load(fp)
    
  # Sanity check of the config file to make sure all the necessary entries exist.
  #
  try:  
    where = "PARAMS"
    cfg['PARAMS']['db_path']
    where = "CONFIGURATIONS"
    for c in cfg['CONFIGURATIONS']:
      c['rotation_data']
      c['error_data']
      c['lens_front_elId']
      c['lens_rear_elId']
      c['mount_front_elId']
      c['mount_rear_elId']
      c['mount_ring_thickness']
      c['error_data_csId']
      c['rotation_data_csId']
      c['hys_idx_1']
      c['hys_idx_2']
      where = "COORDINATE_SYSTEMS"
    for cs in cfg['COORDINATE_SYSTEMS']:
      cs['id']
      cs['angle']
      cs['csMsRecNr']
      where = "DATA"
    for d in cfg['DATA']:
      d['id']
      d['elMsRecNr_range']
      d['elements']
      where = "DATA.elements"
      for el in d['elements']:
	el['id']
	el['elId']
	el['desc']    
  except TypeError:
    print "Configuration file is malformatted, or a key is missing in section " + where

  go(args, cfg)
