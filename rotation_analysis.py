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
  OA_xy_zIs0 = []
  MA_xy_zIs0 = []
  this_elMsRecNr_range = error_data['elMsRecNr_range']
  for elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
    LENS_FRONT_CENTRE_XYZ = []
    LENS_REAR_CENTRE_XYZ = []
    MNT_FRONT_XYZ = []
    MNT_REAR_XYZ = []
    for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
      element = db.getElementFromelRecNr(elRecNr)
      if element['elId'] == configuration['lens_front_elId']:		# front lens
	LENS_FRONT_CENTRE_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
      elif element['elId'] == configuration['lens_rear_elId']:		# rear lens
	LENS_REAR_CENTRE_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
      elif element['elId'] == configuration['mount_front_elId']:	# front mount
	MNT_FRONT_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
      elif element['elId'] == configuration['mount_rear_elId']:		# rear mount
	MNT_REAR_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
 
    if args.oa:
      # optical axis
      OA = axis(LENS_FRONT_CENTRE_XYZ[0], LENS_REAR_CENTRE_XYZ[0])		
      x, y = OA.getXY(z=0)	# evalulate OA at z=0
      OA_xy_zIs0.append((x,y))
    
    if args.ma:
      # mechanical axis
      MA = axis(MNT_FRONT_XYZ[0], MNT_REAR_XYZ[0])
      x, y = MA.getXY(z=0)
      MA_xy_zIs0.append((x,y))

  # Optical axis errors
  #
  if args.oa:
    err = ls_err([xy[0] for xy in OA_xy_zIs0], [xy[1] for xy in OA_xy_zIs0]) 
    OA_r_err_x_y, OA_r_err_euclidean = err.calculate()
  
  # Mechanical axis errors
  #
  if args.ma:
    err = ls_err([xy[0] for xy in MA_xy_zIs0], [xy[1] for xy in MA_xy_zIs0]) 
    MA_r_err_x_y, MA_r_err_euclidean = err.calculate()  

  # Go through each entry in the rotation data and work out the PCS transformed XY position of the 
  # optical and mechanical axes at z=-mount_ring_thickness/2. We also work out the angular deviation of the 
  # axis and the reference (in this case, perpendicular to the mount ring).
  #
  OA_xy_zisMidMountRing = []
  OA_angular_deviation_from_reference = []
  MA_xy_zisMidMountRing = []
  MA_angular_deviation_from_reference = []
  angles = []
  this_elMsRecNr_range = rotation_data['elMsRecNr_range']
  for elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
    LENS_FRONT_CENTRE_XYZ = []
    LENS_REAR_CENTRE_XYZ = []	
    MNT_FRONT_XYZ = []
    MNT_REAR_XYZ = []
    for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
      element = db.getElementFromelRecNr(elRecNr)
      if element['elId'] == configuration['lens_front_elId']:		# front lens
	LENS_FRONT_CENTRE_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
      elif element['elId'] == configuration['lens_rear_elId']:		# rear lens
	LENS_REAR_CENTRE_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
      elif element['elId'] == configuration['mount_front_elId']:	# front mount
	MNT_FRONT_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
      elif element['elId'] == configuration['mount_rear_elId']:		# rear mount
	MNT_REAR_XYZ.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
	
    # optical axis	
    if args.oa:
      OA = axis(LENS_FRONT_CENTRE_XYZ[0], LENS_REAR_CENTRE_XYZ[0])		
      x, y = OA.getXY(z=-(configuration['mount_ring_thickness']/2.))
      angular_deviation = OA.getAngleBetweenOAAndDirectionVector([0,0,1], inDeg=True)
    
      OA_xy_zisMidMountRing.append((x,y))
      OA_angular_deviation_from_reference.append(angular_deviation)
    
    # mechanical axis
    if args.ma:
      MA = axis(MNT_FRONT_XYZ[0], MNT_REAR_XYZ[0])
      x, y = MA.getXY(z=-(configuration['mount_ring_thickness']/2.))
      angular_deviation = MA.getAngleBetweenOAAndDirectionVector([0,0,1], inDeg=True)
    
      MA_xy_zisMidMountRing.append((x,y))
      MA_angular_deviation_from_reference.append(angular_deviation)

    # Get the angle by matching the elMsRecNr with the corresponding entry in the COORDINATE_SYSTEMS 
    # section of the configuration file.
    #
    csys = [c for c in cfg['COORDINATE_SYSTEMS'] if c['csMsRecNr'] == elMsRecNr]
    try:
      assert len(csys) == 1
    except:
      print "Either couldn't find a coordinate system with a csMsRecNr of " + str(elMsRecNr) + ", "\
        "or multiple coordinate systems have been found."
    angles.append(csys[0]['angle'])
  
  # Optical axis sag and hysteresis
  #
  if args.oa:
    sag = ls_sag([xy[0] for xy in OA_xy_zisMidMountRing], [xy[1] for xy in OA_xy_zisMidMountRing])
    OA_r_sag = sag.calculate()
    hys = ls_hys([xy[0] for xy in OA_xy_zisMidMountRing], [xy[1] for xy in OA_xy_zisMidMountRing], angles)
    OA_r_hys = hys.calculate(configuration['hys_idx_1'], configuration['hys_idx_2'])
  
  # Mechanical axis sag and hysteresis
  #
  if args.ma:
    sag = ls_sag([xy[0] for xy in MA_xy_zisMidMountRing], [xy[1] for xy in MA_xy_zisMidMountRing])
    MA_r_sag = sag.calculate()
    hys = ls_hys([xy[0] for xy in MA_xy_zisMidMountRing], [xy[1] for xy in MA_xy_zisMidMountRing], angles)
    MA_r_hys = hys.calculate(configuration['hys_idx_1'], configuration['hys_idx_2'])
  
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
			'angles': (2*np.pi)*np.array(angles)/360.,
			'residuals': OA_r_sag[3]
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
			'angles': (2*np.pi)*np.array(angles)/360.,
			'residuals': MA_r_sag[3]
			}
		      })
    p = p_sag(datasets)
    p.plot()
  
  # And print any information, if requested.
  if args.pi:
    print
    print "[" + configuration['id'] + "]"
    print 
    
    headers = ['AXIS LABEL',
	       'ROTATION ANGLE (degrees)',
	       'XY CENTRE (micron)', 
	       'AXIS DEVIATION FROM REFERENCE (arcmin)'
	       ]
    data = []
    if args.oa:
      for idx, (xy, angle) in enumerate(zip(OA_xy_zisMidMountRing, OA_angular_deviation_from_reference)):
	data.append(['OPTICAL',
		    round(angles[idx]),
		    tuple((round(xy[0]*10**3, 1), round(xy[1]*10**3, 1))), 
		    round(angle*60, 1),
		    ])
    if args.ma:
      for idx, (xy, angle) in enumerate(zip(MA_xy_zisMidMountRing, MA_angular_deviation_from_reference)):
	data.append(['MECHANICAL',
		    round(angles[idx]),
		    tuple((round(xy[0]*10**3, 1), round(xy[1]*10**3, 1))), 
		    round(angle*60, 1),
		    ])
          
    print tabulate(data, headers)
    print 
    
    headers = ['AXIS LABEL',
	       'FIT XY CENTRE (micron)', 
	       'FIT RADIUS (micron)', 
	       'HYSTERESIS (micron)', 
	       'FIT RESIDUAL RMS (micron)'
	       ]
    data = []
    if args.oa:
      data.append(['OPTICAL',
	      tuple((round(OA_r_sag[0]*10**3, 1), round(OA_r_sag[1]*10**3, 1))), 
	      str(round(OA_r_sag[2]*10**3, 1)), 
	      str(round(OA_r_hys*10**3, 1)), 
	      str(round(np.mean(((OA_r_sag[3]*10**3)**2))**0.5,1 ))
	      ])
    if args.ma:
      data.append(['MECHANICAL',
	      tuple((round(MA_r_sag[0]*10**3, 1), round(MA_r_sag[1]*10**3, 1))), 
	      str(round(MA_r_sag[2]*10**3, 1)), 
	      str(round(MA_r_hys*10**3, 1)), 
	      str(round(np.mean(((MA_r_sag[3]*10**3)**2))**0.5,1 ))
	      ])

    print tabulate(data, headers)
    print 

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
      c['csId']
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

      





