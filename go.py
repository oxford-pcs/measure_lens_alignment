#!/usr/local/bin/python

import argparse

import json
import numpy as np
from tabulate import tabulate

from database import CMM_access
from solver import optical_axis_solver
from lens_statistics import sag as ls_sag, hysteresis as ls_hys, measurementError as ls_err
from plotter import sag as p_sag

def solve(args, cfg):
  # Find information regarding this configuration in the config file.
  #
  configuration = [c for c in cfg['CONFIGURATIONS'] if c['id'] == 'lens_1']
  try:
    assert len(configuration) == 1
    configuration = configuration[0]
  except:
    print "Either no configuration, or multiple configurations have been found with the same id."
    exit(0)
  rotation_configuration_data = [d for d in cfg['DATA'] 
				 if d['id'] == configuration['rotation_configuration']]
  error_configuration_data = [d for d in cfg['DATA'] 
			      if d['id'] == configuration['error_configuration']]
  
  # Retrieve the corresponding data for rotation/error calculations.
  #
  try:
    assert len(rotation_configuration_data) == 1
    assert len(error_configuration_data) == 1
    rotation_configuration_data = rotation_configuration_data[0]
    error_configuration_data = error_configuration_data[0]
  except:
    print "Either the datasets specified in the configuration doesn't exist, "\
      "or multiple datasets have been found with the same id."
    exit(0)
    
  try:
    db = CMM_access(cfg['PARAMS']['db_path'])
  except TypeError:
    print "TypeError excepted. Most likely reason is that the database does not exist at the "\
      "location specifed in the configuration file."
   
  # Go through each entry in the error data and work out the PCS transformed XY position of the 
  # optical axis at z=0. We then use this as a measure of the error in both measuring and 
  # fitting the sphere.
  #
  OA_xy_zIs0 = []
  this_elMsRecNr_range = error_configuration_data['elMsRecNr_range']
  for elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
    SPH_1_xy = []
    SPH_2_xy = []
    for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
      element = db.getElementFromelRecNr(elRecNr)
      if element['elId'] == configuration['lens_front_elId']:		# front lens
	SPH_1_xy.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
      elif element['elId'] == configuration['lens_rear_elId']:		# rear lens
	SPH_2_xy.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
    solver = optical_axis_solver(SPH_1_xy[0], SPH_2_xy[0])		
    x, y = solver.getXY(z=0)	# evalulate OA at z=0
    OA_xy_zIs0.append((x,y))

  
  err = ls_err([xy[0] for xy in OA_xy_zIs0], [xy[1] for xy in OA_xy_zIs0]) 
  r_err_x_y, r_err_euclidean = err.calculate()

  # Go through each entry in the rotation data and work out the PCS transformed XY position of the 
  # optical axis at z=-mount_ring_thickness/2. We also work out the angular deviation of the axis 
  # and the reference (in this case, the mount ring).
  #
  OA_xy_zisMidMountRing = []
  OA_angular_deviation_from_reference = []
  angles = []
  this_elMsRecNr_range = rotation_configuration_data['elMsRecNr_range']
  for elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
    SPH_1_xy = []
    SPH_2_xy = []	 
    for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
      element = db.getElementFromelRecNr(elRecNr)
      if element['elId'] == configuration['lens_front_elId']:		# front lens
	SPH_1_xy.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
      elif element['elId'] == configuration['lens_rear_elId']:		# rear lens
	SPH_2_xy.append(db.transElActPos1IntoPCS(element['elRecNr'], configuration['csId']))
    solver = optical_axis_solver(SPH_1_xy[0], SPH_2_xy[0])		
    x, y = solver.getXY(z=0)	# evalulate OA at z=0
    angular_deviation = solver.getAngleBetweenOAAndDirectionVector([0,0,1], inDeg=True)
    
    OA_xy_zisMidMountRing.append((x,y))
    OA_angular_deviation_from_reference.append(angular_deviation)
    
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
  
  sag = ls_sag([xy[0] for xy in OA_xy_zisMidMountRing], [xy[1] for xy in OA_xy_zisMidMountRing])
  r_sag = sag.calculate()

  hys = ls_hys([xy[0] for xy in OA_xy_zisMidMountRing], [xy[1] for xy in OA_xy_zisMidMountRing], angles)
  r_hys = hys.calculate(configuration['hysIdx1'], configuration['hysIdx2'])
  
  # Now we can plot, if requested. We construct datasets first in case we want to plot optical and 
  # mechanical results on the same axes.
  #
  if args.p2d:
    datasets = []
    datasets.append({'heading': 'Optical Axis',
		     'data': {
		       'x': [xy[0] for xy in OA_xy_zisMidMountRing],
		       'y': [xy[1] for xy in OA_xy_zisMidMountRing],
		       'x_err': np.amax(r_err_x_y),
		       'y_err': np.amax(r_err_x_y),
		       'fit_xc': r_sag[0],
		       'fit_yc': r_sag[1],
		       'fit_r': r_sag[2],
		       'angles': (2*np.pi)*np.array(angles)/360.,
		       'residuals': r_sag[3]
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
    for idx, (xy, angle) in enumerate(zip(OA_xy_zisMidMountRing, OA_angular_deviation_from_reference)):
      data.append(['OPTICAL',
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
    data = [['OPTICAL',
	     tuple((round(r_sag[0]*10**3, 1), round(r_sag[1]*10**3, 1))), 
	     str(round(r_sag[2]*10**3, 1)), 
	     str(round(r_hys*10**3, 1)), 
	     str(round(np.mean(((r_sag[3]*10**3)**2))**0.5,1 ))
	     ]]

    print tabulate(data, headers)
    print 

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-l", help="configuration id to process", default="lens_1")
  parser.add_argument("-pi", help="print info?", action='store_true')
  parser.add_argument("-p2d", help="plot 2d?", action='store_true')
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
      c['rotation_configuration']
      c['error_configuration']
      c['lens_front_elId']
      c['lens_rear_elId']
      c['mount_lens_thickness']
      c['csId']
      c['hysIdx1']
      c['hysIdx2']
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

  solve(args, cfg)

      





