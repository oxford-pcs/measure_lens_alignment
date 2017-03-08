#!/usr/local/bin/python

import argparse

import json
import numpy as np
from tabulate import tabulate

from database import CMM_access
from solver import axis
from lens_statistics import sag as ls_sag, hysteresis as ls_hys, measurementError as ls_err
from plotter import sag as p_sag

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
   
  # Go through each entry in the error data (multiple measurements of single rotation positions) and work 
  # out the PCS transformed XY position of the optical and mechanical axes at z=0. We then use this as a 
  # measure of the error in both measuring and fitting.
  #
  OA_err_axes = []
  MA_err_axes = []
  this_elMsRecNr_range = error_data['elMsRecNr_range']
  for elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
    LENS_FRONT_CENTRE_XYZ = []
    LENS_REAR_CENTRE_XYZ = []
    MNT_FRONT_XYZ = []
    MNT_REAR_XYZ = []
    for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
      element = db.getElementFromelRecNr(elRecNr)
      if element['elId'] == configuration['lens_front_elId']:		# front lens
	LENS_FRONT_CENTRE_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['error_data_csId']))
      elif element['elId'] == configuration['lens_rear_elId']:		# rear lens
	LENS_REAR_CENTRE_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['error_data_csId']))
      elif element['elId'] == configuration['mount_front_elId']:	# front mount
	MNT_FRONT_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['error_data_csId']))
      elif element['elId'] == configuration['mount_rear_elId']:		# rear mount
	MNT_REAR_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['error_data_csId']))
 
    if args.oa:
      # optical axis
      try:
        OA_err_axes.append(axis(LENS_FRONT_CENTRE_XYZ[0], LENS_REAR_CENTRE_XYZ[0]))
      except IndexError:
	print "Optical axis error data is empty."
	exit(0)

    if args.ma:
      # mechanical axis
      try:
        MA_err_axes.append(axis(MNT_FRONT_XYZ[0], MNT_REAR_XYZ[0]))
      except IndexError:
	print "Mechanical axis error data is empty."
	exit(0)


  # Optical axis errors
  #
  if args.oa:
    err = ls_err(OA_err_axes) 
    OA_r_err_x_y, OA_r_err_euclidean = err.calculate_position_error_at_z(z=-(configuration['mount_ring_thickness']/2.))
    OA_err_angle = err.calculate_angle_error_at_z(z=-(configuration['mount_ring_thickness']/2.))
    
  # Mechanical axis errors
  #
  if args.ma:
    err = ls_err(MA_err_axes) 
    MA_r_err_x_y, MA_r_err_euclidean = err.calculate_position_error_at_z(z=-(configuration['mount_ring_thickness']/2.)) 
    MA_err_angle = err.calculate_angle_error_at_z(z=-(configuration['mount_ring_thickness']/2.))
    
  # Go through each entry in the rotation data and work out the PCS transformed XY position of the 
  # optical and mechanical axes at z=-mount_ring_thickness/2. We also work out the angular deviation of the 
  # axis and the reference (in this case, perpendicular to the mount ring).
  #
  optical_axes = []
  mechanical_axes = []
  mount_angles = []
  this_elMsRecNr_range = rotation_data['elMsRecNr_range']
  for elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
    LENS_FRONT_CENTRE_XYZ = []
    LENS_REAR_CENTRE_XYZ = []	
    MNT_FRONT_XYZ = []
    MNT_REAR_XYZ = []
    for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
      element = db.getElementFromelRecNr(elRecNr)
      if element['elId'] == configuration['lens_front_elId']:		# front lens
	LENS_FRONT_CENTRE_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['rotation_data_csId']))
      elif element['elId'] == configuration['lens_rear_elId']:		# rear lens
	LENS_REAR_CENTRE_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['rotation_data_csId']))
      elif element['elId'] == configuration['mount_front_elId']:	# front mount
	MNT_FRONT_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['rotation_data_csId']))
      elif element['elId'] == configuration['mount_rear_elId']:		# rear mount
	MNT_REAR_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['rotation_data_csId']))
     
    # optical axis	
    if args.oa:
      try:
        OA = axis(LENS_FRONT_CENTRE_XYZ[0], LENS_REAR_CENTRE_XYZ[0])		
      except IndexError:
	print "Optical axis error data is empty."
	exit(0)
      optical_axes.append(OA)
    
    # mechanical axis
    if args.ma:
      try:
        MA = axis(MNT_FRONT_XYZ[0], MNT_REAR_XYZ[0])
      except IndexError:
	print "Optical axis error data is empty."
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
    
  # Optical axis sag, angular deviation from [0, 0, 1] and hysteresis
  #
  OA_xy_zisMidMountRing = []
  OA_angles_from_mount_normal = []
  OA_xy_angles_from_12_o_clock = []
  if args.oa:
    # get XY at z=0 and the angle between the optical axis and the mount vector [0, 0, 1] (i.e. z axis)
    for ax in optical_axes:
      OA_xy_zisMidMountRing.append((ax.getXY(z=-(configuration['mount_ring_thickness']/2.))))
      OA_angles_from_mount_normal.append(ax.getComponentAngles(inDeg=True))
      
    # sag  
    sag = ls_sag([xy[0] for xy in OA_xy_zisMidMountRing], [xy[1] for xy in OA_xy_zisMidMountRing])
    OA_r_sag = sag.calculate()
    
    # hysteresis
    hys = ls_hys([xy[0] for xy in OA_xy_zisMidMountRing], [xy[1] for xy in OA_xy_zisMidMountRing], mount_angles)
    OA_r_hys = hys.calculate(configuration['hys_idx_1'], configuration['hys_idx_2'])
    
    # determine the angular rotation of the optical axis xy vector relative to the reference vector [0, 1] 
    # (clockwise from 12 o'clock) using the predetermined fit centre.
    for ax in optical_axes:  
      ref_axis = np.array([0, 1])	# unit vector reference axis for angle
      point_vector_from_fit_centre = np.array((ax.getXY(z=-(configuration['mount_ring_thickness']/2.)))) - OA_r_sag[0:2]
    
      dotP = np.dot(ref_axis, point_vector_from_fit_centre)
      crossP = np.cross(ref_axis, point_vector_from_fit_centre)
      angle = np.arccos(dotP/(np.linalg.norm(ref_axis)*np.linalg.norm(point_vector_from_fit_centre)))
      if np.sign(crossP) > 0:
	angle = (np.pi-angle) + np.pi
      OA_xy_angles_from_12_o_clock.append((360*angle)/(2*np.pi))
  
  # Mechanical axis sag, angular deviation from [0, 0, 1] and hysteresis
  #
  MA_xy_zisMidMountRing = []
  MA_angles_from_mount_normal = []
  MA_xy_angles_from_12_o_clock = []
  if args.ma:
    # get XY at z=0 and the angle between the mechanical axis and the mount vector [0, 0, 1] (i.e. z axis)
    for ax in mechanical_axes:
      MA_xy_zisMidMountRing.append((ax.getXY(z=-(configuration['mount_ring_thickness']/2.))))
      MA_angles_from_mount_normal.append(ax.getComponentAngles(inDeg=True))
      
    # sag  
    sag = ls_sag([xy[0] for xy in MA_xy_zisMidMountRing], [xy[1] for xy in MA_xy_zisMidMountRing])
    MA_r_sag = sag.calculate()
    
    # hysteresis
    hys = ls_hys([xy[0] for xy in MA_xy_zisMidMountRing], [xy[1] for xy in MA_xy_zisMidMountRing], mount_angles)
    MA_r_hys = hys.calculate(configuration['hys_idx_1'], configuration['hys_idx_2'])
    
    # determine the angular rotation of the mechanical axis xy vector relative to the reference vector [0, 1]
    # (clockwise from 12 o'clock) using the predetermined fit centre.
    for ax in mechanical_axes:  
      ref_axis = np.array([0, 1])	# unit vector reference axis for angle
      point_vector_from_fit_centre = np.array((ax.getXY(z=-(configuration['mount_ring_thickness']/2.)))) - MA_r_sag[0:2]
    
      dotP = np.dot(ref_axis, point_vector_from_fit_centre)
      crossP = np.cross(ref_axis, point_vector_from_fit_centre)
      angle = np.arccos(dotP/(np.linalg.norm(ref_axis)*np.linalg.norm(point_vector_from_fit_centre)))
      if np.sign(crossP) > 0:
	angle = (np.pi-angle) + np.pi
      MA_xy_angles_from_12_o_clock.append((360*angle)/(2*np.pi))
    
  # Print any information, if requested.
  if args.pi:
    print
    print "[" + configuration['id'] + "]"
    print 
    
    headers = ['AXIS LABEL',
	       'MOUNT POSITION',
	       'XY CENTRE (micron)', 
	       'AXIS ANGLE FROM Z AXIS (arcmin)',
	       'XY ANGLE FROM 12 o\'clock'
	       ]
    
    data = []
    if args.oa:
      for (mount_angle, xy, mount_normal_angle, xy_angle) in zip(mount_angles, OA_xy_zisMidMountRing, OA_angles_from_mount_normal, OA_xy_angles_from_12_o_clock):
	data.append(['OPTICAL',
		    round(mount_angle),
		    tuple((round(xy[0]*10**3, 1), round(xy[1]*10**3, 1))), 
		    round(mount_normal_angle[2]*60, 2),
		    round(xy_angle)
		    ])
    if args.ma:
      for (mount_angle, xy, mount_normal_angle, xy_angle) in zip(mount_angles, MA_xy_zisMidMountRing, MA_angles_from_mount_normal, MA_xy_angles_from_12_o_clock):
	data.append(['MECHANICAL',
		    round(mount_angle),
		    tuple((round(xy[0]*10**3, 1), round(xy[1]*10**3, 1))), 
		    round(mount_normal_angle[2]*60, 2),
		    round(xy_angle)
		    ])
          
    print tabulate(data, headers)
    print 
    
    headers = ['AXIS LABEL',
	       'ERROR POS X (micron)',
	       'ERROR POS Y (micron)',	  
	       'FIT XY CENTRE (micron)', 
	       'FIT RADIUS (micron)', 
	       'HYSTERESIS (micron)', 
	       'FIT RESIDUAL RMS (micron)',    
	       'ERROR Z AXIS ANGLE (arcmin)',
	       ]
    data = []
    if args.oa:
      data.append(['OPTICAL',
	      round(OA_r_err_x_y[0]*10**3, 1),
	      round(OA_r_err_x_y[1]*10**3, 1),
	      tuple((round(OA_r_sag[0]*10**3, 1), round(OA_r_sag[1]*10**3, 1))), 
	      str(round(OA_r_sag[2]*10**3, 1)), 
	      str(round(OA_r_hys*10**3, 1)), 
	      str(round(np.mean(((OA_r_sag[3]*10**3)**2))**0.5,1)),
	      round(OA_err_angle*60, 2)
	      ])
    if args.ma:
      data.append(['MECHANICAL',
	      round(MA_r_err_x_y[0]*10**3, 1),
	      round(MA_r_err_x_y[1]*10**3, 1),
	      tuple((round(MA_r_sag[0]*10**3, 1), round(MA_r_sag[1]*10**3, 1))), 
	      str(round(MA_r_sag[2]*10**3, 1)), 
	      str(round(MA_r_hys*10**3, 1)), 
	      str(round(np.mean(((MA_r_sag[3]*10**3)**2))**0.5,1)),
	      round(MA_err_angle*60, 2)
	      ])

    print tabulate(data, headers)
    print     
  
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
			'x_err': np.amax(OA_r_err_x_y),
			'y_err': np.amax(OA_r_err_x_y),
			'fit_xc': OA_r_sag[0],
			'fit_yc': OA_r_sag[1],
			'fit_r': OA_r_sag[2],
			'mount_angles': (2*np.pi)*np.array(mount_angles)/360.,
			'residuals': OA_r_sag[3],
			'angles_from_mount_normal': OA_angles_from_mount_normal,
			'xy_angles_from_12_o_clock': (2*np.pi)*np.array(OA_xy_angles_from_12_o_clock)/360.
			}
		      })
    if args.ma:
      datasets.append({'heading': 'Mechanical Axis',
		      'data': {
			'x': [xy[0] for xy in MA_xy_zisMidMountRing],
			'y': [xy[1] for xy in MA_xy_zisMidMountRing],
			'x_err': np.amax(MA_r_err_x_y),
			'y_err': np.amax(MA_r_err_x_y),
			'fit_xc': MA_r_sag[0],
			'fit_yc': MA_r_sag[1],
			'fit_r': MA_r_sag[2],
			'mount_angles': (2*np.pi)*np.array(mount_angles)/360.,
			'residuals': MA_r_sag[3],
			'angles_from_mount_normal': MA_angles_from_mount_normal,
			'xy_angles_from_12_o_clock': (2*np.pi)*np.array(MA_xy_angles_from_12_o_clock)/360.
			}
		      })
    p = p_sag(datasets, hard=True)
    p.plot()
 
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-l", help="configuration id to process", default="lens_1")
  parser.add_argument("-pi", help="print info?", action='store_true')
  parser.add_argument("-p2d", help="plot 2d?", action='store_true')
  parser.add_argument("-oa", help="consider optical axis?", action='store_true')
  parser.add_argument("-ma", help="consider mechanical axis?", action='store_true')
  args = parser.parse_args()
 
  CONFIG_FILE = "config.json"
  with open("config.json") as fp:
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

      





