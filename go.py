#!/usr/local/bin/python

import pymdb
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from tabulate import tabulate
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch, Rectangle
from mpl_toolkits.mplot3d import proj3d
import argparse
from scipy.optimize import curve_fit

class Arrow3D(FancyArrowPatch):
  def __init__(self, xs, ys, zs, *args, **kwargs):
    FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
    self._verts3d = xs, ys, zs
  def draw(self, renderer):
    xs3d, ys3d, zs3d = self._verts3d
    xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
    self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
    FancyArrowPatch.draw(self, renderer)

def leastsq_circle(x,y):
  def calc_R(x,y, xc, yc):
      """ calculate the distance of each 2D points from the center (xc, yc) """
      return np.sqrt((x-xc)**2 + (y-yc)**2)

  def f(c, x, y):
      """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
      Ri = calc_R(x, y, *c)
      return Ri - Ri.mean()

  # coordinates of the barycenter
  x_m = np.mean(x)
  y_m = np.mean(y)
  center_estimate = x_m, y_m
  center, ier = optimize.leastsq(f, center_estimate, args=(x,y))
  xc, yc = center
  Ri       = calc_R(x, y, *center)
  R        = Ri.mean()
  residu   = (Ri - R)**2
  return xc, yc, R, residu
	    
def gaus(x,a,x0,sigma):
  return a*np.exp(-(x-x0)**2/(2*sigma**2))

def solve(lens_configuration_name, args, DATA_ONLY=False, opt_f_error=None, opt_r_error=None):
  MECH_F = []		# front mechanical lens diameter
  MECH_R = []		# rear mechanical lens diameter
  OPT_F = []		# front lens surface
  OPT_R = []		# rear lens surface
  MECH = []		# front mechanical lens cylinder
  OPT = []		# optical axis

  print [d for d in cfg['DATA']]
  exit(0)
  this_lens_data = None
  for d in cfg['DATA']:			# which lens configuration are we processing?
      if d['id'] == lens_configuration_name:
	  this_lens_data = d
  if this_lens_data == None:      	# we couldn't find the configuration
      return None

  this_elMsRecNr_range = this_lens_data['elMsRecNr_range']									# get the range of elMsRecNr that correspond to this set of measurements
  for entry in db[t_el]['data']:												# for each element in the elements table
      this_elId = entry[db[t_el]['headers'].index('elId')]									# get the element ID (e.g. CIR_1)
      this_elMsRecNr = int(entry[db[t_el]['headers'].index('elMsRecNr')])							# get the element measurement number
      for el in this_lens_data['elements']:											# scan through all the elements specified in the config file
	  if el['elId'] == this_elId and this_elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):		# match it, provided it's within the range
	      for csys in cfg['COORDINATE_SYSTEMS']:										# find the corresponding coordinate system			
		  if csys['csMsRecNr'] == this_elMsRecNr:
		      try:
			  t_inv = csys['t_inv']
			  x1 = float(entry[db[t_el]['headers'].index('elActPos1X')]) 		# elActPos1* are evaluated at z=0
			  y1 = float(entry[db[t_el]['headers'].index('elActPos1Y')])
			  z1 = float(entry[db[t_el]['headers'].index('elActPos1Z')])
			  x1y1z1_column_vector = np.array([[x1], [y1], [z1], [1]])	
			  x1y1z1_transformed = np.dot(t_inv, x1y1z1_column_vector)		# PCS coordinates 

			  elActDim1 = float(entry[db[t_el]['headers'].index('elActDim1')])	# "Dim" is radius in the case of a sphere or circle
			  
			  if el['key'] == 'MECH_F' or el['key'] == 'MECH_R' or \
			      el['key'] == 'OPT_F' or el['key'] == 'OPT_R':
			      x1t = round(x1y1z1_transformed[0]*10**3, 3)	# micron
			      y1t = round(x1y1z1_transformed[1]*10**3, 3)
			      z1t = round(x1y1z1_transformed[2]*10**3, 3) 
			      r = np.sqrt((x1t**2) + (y1t**2) + (z1t**2))
			      theta = 2*np.pi*float(csys['angle'])/360
			      if el['key'] == 'MECH_F':
				  MECH_F.append({'x': x1t, 'y': y1t, 'z': z1t, 'theta': theta, 'r': r, 'radius': elActDim1})
			      elif el['key'] == 'MECH_R':
				  MECH_R.append({'x': x1t, 'y': y1t, 'z': z1t, 'theta': theta, 'r': r, 'radius': elActDim1}) 
			      elif el['key'] == 'OPT_F':
				  OPT_F.append({'x': x1t, 'y': y1t, 'z': z1t, 'theta': theta, 'r': r, 'radius': elActDim1})
			      elif el['key'] == 'OPT_R':
				  OPT_R.append({'x': x1t, 'y': y1t, 'z': z1t, 'theta': theta, 'r': r, 'radius': elActDim1})
		      except KeyError:
			  exit(0)	

  # find optical axis at z=0
  '''
      x, y, z = xyz_f(x, y, z) + t*dir_v(x, y, z)

      where dir_v = xyz_f - xyz_r
      
      parametric form:
      
      z = xyz_f.z + t*dir_v.z

      solve for z=0:
      
      t = -xyz_f.z/dir_v.z
  '''	
  idx_angle = 0
  for f, r in zip(OPT_F, OPT_R):
      if f['theta'] != r['theta']:								# should always be in same angle order, but check anyway
	  exit(0)
      xyz_f = np.array([f['x'], f['y'], f['z']])
      xyz_r = np.array([r['x'], r['y'], r['z']])
      if opt_f_error is not None:
	xyz_f = xyz_f + opt_f_error[idx_angle]
      if opt_r_error is not None:
	xyz_r = xyz_r + opt_r_error[idx_angle]

      dir_v = xyz_f - xyz_r									# direction vector of OA
      dir_v_n = dir_v/np.linalg.norm(dir_v)							# normalise OA direction vector
      zis0 = -xyz_f[2]/dir_v[2]
      oa_at_zis0 = xyz_f + zis0*dir_v
      x = oa_at_zis0[0]
      y = oa_at_zis0[1]
      z = oa_at_zis0[2]
      R = np.sqrt((x**2) + (y**2) + (z**2))
      
      # direction cosines
      cosX = dir_v_n[0]/np.sqrt((dir_v_n[0]**2)+(dir_v_n[1]**2)+(dir_v_n[2]**2))
      cosY = dir_v_n[1]/np.sqrt((dir_v_n[0]**2)+(dir_v_n[1]**2)+(dir_v_n[2]**2))
      cosZ = dir_v_n[2]/np.sqrt((dir_v_n[0]**2)+(dir_v_n[1]**2)+(dir_v_n[2]**2))
      x_angle = np.arccos(cosX)									# component angle between dir_v and x axis
      y_angle = np.arccos(cosY)
      z_angle = np.arccos(cosZ)  
      
      xyz_lens_surface_axis = np.array([0, 0, f['z']])						# normal to lens surface
      xyz_lens_surface_axis_n = xyz_lens_surface_axis						# normalise lens surface axis direction vector

      dotP = np.dot(dir_v, xyz_lens_surface_axis)
	
      OA_angle_lens_surface = np.arccos(dotP/(np.linalg.norm(dir_v)*np.linalg.norm(xyz_lens_surface_axis))) # find angle from dot product = |dir_v_n|*|xyz_lens_surface_axis_n|*cos(angle)

      OPT.append({'OA_reference_angle_x': np.abs((np.pi/2)-x_angle), 				# component angles
		  'OA_reference_angle_y': np.abs((np.pi/2)-y_angle),
		  'OA_reference_angle_z': np.abs((np.pi)-z_angle),					
		  'x': x, 'y': y, 'z': z, 'theta': f['theta'], 					# at z=0
		  'r': R, 'OA_reference_angle': OA_angle_lens_surface
	    })  
      idx_angle = idx_angle+1

  # CALCULATIONS
  # ------------
  # mech_f
  x, y = [i['x'] for i in MECH_F], [i['y'] for i in MECH_F]
  mech_f_xc, mech_f_yc, mech_f_R, mech_f_residu = leastsq_circle([i['x'] for i in MECH_F], [i['y'] for i in MECH_F])
  theta_fit = np.linspace(-np.pi, np.pi, 180)
  mech_f_x_fit = mech_f_xc + mech_f_R*np.cos(theta_fit)
  mech_f_y_fit = mech_f_yc + mech_f_R*np.sin(theta_fit)    
  mech_f_xd = x[-1]-x[0]
  mech_f_yd = y[-1]-y[0]
  mech_f_d = np.sqrt((mech_f_xd**2)+(mech_f_yd**2))    
  
  # mech_r 
  x, y = [i['x'] for i in MECH_R], [i['y'] for i in MECH_R]
  mech_r_xc, mech_r_yc, mech_r_R, mech_r_residu = leastsq_circle([i['x'] for i in MECH_R], [i['y'] for i in MECH_R])
  theta_fit = np.linspace(-np.pi, np.pi, 180)
  mech_r_x_fit = mech_r_xc + mech_r_R*np.cos(theta_fit)
  mech_r_y_fit = mech_r_yc + mech_r_R*np.sin(theta_fit)
  mech_r_xd = x[-1]-x[0]
  mech_r_yd = y[-1]-y[0]
  mech_r_d = np.sqrt((mech_r_xd**2)+(mech_r_yd**2))
  
  # opt
  x, y = [i['x'] for i in OPT], [i['y'] for i in OPT]
  opt_xc, opt_yc, opt_R, opt_residu = leastsq_circle([i['x'] for i in OPT], [i['y'] for i in OPT])
  theta_fit = np.linspace(-np.pi, np.pi, 180)
  opt_x_fit = opt_xc + opt_R*np.cos(theta_fit)
  opt_y_fit = opt_yc + opt_R*np.sin(theta_fit)
  opt_xd = x[-1]-x[0]
  opt_yd = y[-1]-y[0]
  opt_d = np.sqrt((opt_xd**2)+(opt_yd**2))

  if DATA_ONLY:
    return MECH_F, MECH_R, OPT_F, OPT_R, OPT
  
  if args.p2d:
    nplots = 4
    f, axes = plt.subplots(nplots,1, figsize=(22,9))
    '''
      Produce axis plots (mech and opt) in both linear and polar coordinates
    '''
    ax = plt.subplot(1, nplots, 1)
    
    # LINEAR DISPLACEMENTS
    ## MECH_F
    x, y = [i['x'] for i in MECH_F], [i['y'] for i in MECH_F]
    ax.plot(x, y, linewidth=3, label='MECHANICAL FRONT', color='r', ls='-')
    for i in MECH_F:
	ax.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    ax.plot(mech_f_x_fit, mech_f_y_fit, color='r', ls='--', label='Circular fit')
    ax.plot(mech_f_xc, mech_f_yc, color='r', marker='x')

    ## MECH_R
    x, y = [i['x'] for i in MECH_R], [i['y'] for i in MECH_R]
    ax.plot(x, y, linewidth=3, label='MECHANICAL REAR', color='b', ls='-') 
    for i in MECH_R:
	ax.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    ax.plot(mech_r_x_fit, mech_r_y_fit, color='b', ls='--', label='Circular fit')
    ax.plot(mech_r_xc, mech_r_yc, color='g', marker='x')

    ## OPT
    x, y = [i['x'] for i in OPT], [i['y'] for i in OPT]
    ax.plot(x, y, linewidth=3, label='OPTICAL AXIS', color='g', ls='-')
    for i in OPT:
	ax.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    ax.plot(opt_x_fit, opt_y_fit, color='g', ls='--', label='Circular fit')
    ax.plot(opt_xc, opt_yc, color='b', marker='x')

    ax.legend(loc='lower right', prop={'size':8})
    ax.set_title("Decentre relative to lens ring mount")
    ax.set_xlabel("x (micron)")
    ax.set_ylabel("y (micron)")
    ax.set(aspect=1, adjustable='datalim')

    # RADIAL DISPLACEMENTS
    ax = plt.subplot(1, nplots, 2, projection='polar')
    ax.plot([i['theta'] for i in MECH_F], [i['r'] for i in MECH_F], linewidth=3, color='r', ls='-', label='MECHANICAL FRONT')
    ax.plot([i['theta'] for i in MECH_R], [i['r'] for i in MECH_R], linewidth=3, color='g', ls='-', label='MECHANICAL REAR')
    ax.plot([i['theta'] for i in OPT], [i['r'] for i in OPT], linewidth=3, color='b', ls='-', label='OPTICAL AXIS')
    (lines,labels) = plt.thetagrids( range(0,360,60), range(0,360,60))

    # RESIDUALS FROM CIRCULAR FIT
    ax = plt.subplot(1, nplots, 3, projection='polar')
    ax.plot([i['theta'] for i in MECH_F], mech_f_residu, color='r', ls='-', label='MECHANICAL FRONT')
    ax.plot([i['theta'] for i in MECH_R], mech_r_residu, color='g', ls='-', label='MECHANICAL REAR')
    ax.plot([i['theta'] for i in OPT], opt_residu, color='b', ls='-', label='OPTICAL AXIS')
    ax.set_title("Residual from circular fit (micron)")
    (lines,labels) = plt.thetagrids( range(0,360,60), range(0,360,60))
	    
    # OPTICAL AXIS REFERENCE ANGLE RECESSION
    ax = plt.subplot(1, nplots, 4, projection='polar')
    ax.plot([i['theta'] for i in OPT], [(360*i['OA_reference_angle'])/(2*np.pi) for i in OPT], color='b', ls='-', label='OPTICAL AXIS')
    ax.plot([i['theta'] for i in OPT], [(360*i['OA_reference_angle_x'])/(2*np.pi) for i in OPT], color='b', lw=2, ls='--', label='OPTICAL AXIS (x)')
    ax.plot([i['theta'] for i in OPT], [(360*i['OA_reference_angle_y'])/(2*np.pi) for i in OPT], color='b', lw=2, ls=':', label='OPTICAL AXIS (y)')
    plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4), ncol=1, prop={'size':12})
    ax.set_title("Angular deviation from coordinate system axis (deg)")
    (lines,labels) = plt.thetagrids( range(0,360,60), range(0,360,60))
    
    plt.show()
    
  if args.p3d:
    f_c = []	# front centres (mm)
    r_c = []
    f_c_interped = []
    r_c_interped = []
    OA_angle_lens_surface_interped = []
    angles = []
    angles_interped = []
    for this_of, this_or, this_oa in zip(OPT_F, OPT_R, OPT):
      f_c.append(np.array([this_of['x']/1000,this_of['y']/1000,this_of['z']/1000]))
      r_c.append(np.array([this_or['x']/1000,this_or['y']/1000,this_or['z']/1000]))
      angles.append(this_oa['theta']*360/(2*np.pi))
      
    max_f_c_z = np.mean([xyz[2] for xyz in f_c])
    max_r_c_z = np.mean([xyz[2] for xyz in r_c])
    f_c_z_normalised = []
    for entry in f_c:
      f_c_z_normalised.append(np.array([entry[0], entry[1], entry[2]-max_f_c_z]))
    r_c_z_normalised = []
    for entry in r_c:
      r_c_z_normalised.append(np.array([entry[0], entry[1], 1+(entry[2]-max_r_c_z)]))
    
    n_interp_steps = 60
    angle_idx = 0
    min_x_f, max_x_f = np.min([xyz[0] for xyz in f_c]), np.max([xyz[0] for xyz in f_c])	# need consistent axes
    min_y_f, max_y_f = np.min([xyz[1] for xyz in f_c]), np.max([xyz[1] for xyz in f_c])
    min_z_f, max_z_f = np.min([xyz[2] for xyz in f_c_z_normalised]), np.max([xyz[2] for xyz in f_c_z_normalised])
    min_x_r, max_x_r = np.min([xyz[0] for xyz in r_c]), np.max([xyz[0] for xyz in r_c])
    min_y_r, max_y_r = np.min([xyz[1] for xyz in r_c]), np.max([xyz[1] for xyz in r_c])
    min_z_r, max_z_r = np.min([xyz[2] for xyz in r_c_z_normalised]), np.max([xyz[2] for xyz in r_c_z_normalised])
    max_OA_reference_angle = np.max(np.array([i['OA_reference_angle'] for i in OPT])*360/(2*np.pi))
    
    fnames = []
    for ii, jj, kk, mm in zip(f_c[:-1], f_c[1:], r_c[:-1], r_c[1:]):
      ii_z_normalised = np.array([ii[0], ii[1], ii[2]-max_f_c_z])
      jj_z_normalised = np.array([jj[0], jj[1], jj[2]-max_f_c_z])
      kk_z_normalised = np.array([kk[0], kk[1], 1+(kk[2]-max_r_c_z)])
      mm_z_normalised = np.array([mm[0], mm[1], 1+(mm[2]-max_r_c_z)])
      increment_vector_f = (jj_z_normalised-ii_z_normalised)/n_interp_steps
      increment_vector_r = (mm_z_normalised-kk_z_normalised)/n_interp_steps
      for i in range(0, n_interp_steps):
	xyz_f = np.array(ii+(i*increment_vector_f))
        xyz_r = np.array(kk+(i*increment_vector_r))
        xyz_f_z_normalised = np.array(ii_z_normalised+(i*increment_vector_f))
        xyz_r_z_normalised = np.array(kk_z_normalised+(i*increment_vector_r))
        f_c_interped.append(xyz_f_z_normalised)
        r_c_interped.append(xyz_r_z_normalised)
        angular_increment = (angles[0]+angles[1])/n_interp_steps
        this_interped_angle = angles[angle_idx]+(angular_increment*i)
        angles_interped.append(this_interped_angle)
        
        fig = plt.figure(figsize=plt.figaspect(2))
        ax = fig.add_subplot(2, 1, 1, projection='3d')
        xlim = [np.min([min_x_f, min_x_r]), np.max([max_x_f, max_x_r])]
        ylim = [np.min([min_y_f, min_y_r]), np.max([max_y_f, max_y_r])]
        zlim = [np.min([min_z_f, min_z_r]), np.max([max_z_f, max_z_r])]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_zlim(zlim)
        ax.view_init(elev=30., azim=45)

        ax.scatter([xyz_f_z_normalised[0]],[xyz_f_z_normalised[1]],[xyz_f_z_normalised[2]], color="r", s=10)
        ax.scatter([xyz_r_z_normalised[0]],[xyz_r_z_normalised[1]],[xyz_r_z_normalised[2]], color="b", s=10)

        arrow_OA = Arrow3D([xyz_f_z_normalised[0],xyz_r_z_normalised[0]],[xyz_f_z_normalised[1],xyz_r_z_normalised[1]],
			   [xyz_f_z_normalised[2],xyz_r_z_normalised[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color='k')
        ax.add_artist(arrow_OA)
        
        arrow_REFCOORDSYS = Arrow3D([0,0],[0,0],zlim, mutation_scale=20, lw=1, arrowstyle="-|>", color='k', ls='--')
        ax.add_artist(arrow_REFCOORDSYS)
        
        ax.plot([xyz[0] for xyz in f_c_interped],[xyz[1] for xyz in f_c_interped], [zlim[0] for xyz in f_c_interped], 'r-')
        ax.plot([xyz[0] for xyz in f_c_interped],[ylim[0] for xyz in f_c_interped], [xyz[2] for xyz in f_c_interped], 'r-')
        ax.plot([xlim[0] for xyz in f_c_interped], [xyz[1] for xyz in f_c_interped], [xyz[2] for xyz in f_c_interped], 'r-')
        ax.plot([xyz[0] for xyz in r_c_interped],[xyz[1] for xyz in r_c_interped], [zlim[0] for xyz in r_c_interped], 'b-')
        ax.plot([xyz[0] for xyz in r_c_interped],[ylim[0] for xyz in r_c_interped], [xyz[2] for xyz in r_c_interped], 'b-')
        ax.plot([xlim[0] for xyz in r_c_interped], [xyz[1] for xyz in r_c_interped], [xyz[2] for xyz in r_c_interped], 'b-')
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        ax.annotate(str(int(round(this_interped_angle,0))) + ' deg', xy=(0.7,0.95), xycoords='figure fraction')
        
        #FIXME: this should be a function
        dir_v = xyz_f - xyz_r				# direction vector of OA
	dir_v_n = dir_v/np.linalg.norm(dir_v)		# normalise OA direction vector
	zis0 = -xyz_f[2]/dir_v[2]
	oa_at_zis0 = xyz_f + zis0*dir_v
	x = oa_at_zis0[0]
	y = oa_at_zis0[1]
	z = oa_at_zis0[2]
	R = np.sqrt((x**2) + (y**2) + (z**2))
	
	xyz_lens_surface_axis = np.array([0, 0, zlim[0]])					# normal to lens surface
        xyz_lens_surface_axis_n = xyz_lens_surface_axis/np.linalg.norm(xyz_lens_surface_axis)	# normalise lens surface axis direction vector
	
	dotP = np.dot(dir_v, xyz_lens_surface_axis)
	
        OA_angle_lens_surface = np.arccos(dotP/(np.linalg.norm(dir_v)*np.linalg.norm(xyz_lens_surface_axis))) # find angle from dot product = |dir_v_n|*|xyz_lens_surface_axis_n|*cos(angle)

	OA_angle_lens_surface_interped.append((OA_angle_lens_surface*360/(2*np.pi)))
	
        ax = fig.add_subplot(3, 1, 3, projection='polar')
        ax.plot(np.array(angles_interped)*(2*np.pi)/360, OA_angle_lens_surface_interped, linewidth=3, color='b', ls='-', label='OPTICAL AXIS')
        ax.set_rmax(max_OA_reference_angle)
        (lines,labels) = plt.thetagrids( range(0,360,60), range(0,360,60))
        
        fname = str(int(round(angles[angle_idx],2))) + '_' + str(i) + ".png"
        fnames.append(fname)
        plt.savefig(fname)
        plt.close()  

      angle_idx = angle_idx + 1  
    
    script_text = '#!/bin/bash\n/usr/bin/convert -delay 5 ' + ' '.join(fnames) + ' -loop 0 animated.gif'
    with open("makegif.csh", 'w') as f:
      f.write(script_text)
    
  if args.pi:
    '''
      Pretty print pertinent data
    '''
    
    print OPT
        
    headers = ['Axis', 'AXIS MAX DEVIATION (deg)', 'FIT_XCENT (micron)', 'FIT_YCENT (micron)', 'FIT_RADIUS (micron)', 'ELASTICITY (micron)', 'FIT_RESIDUAL_RMS (micron)']
    data = [['MECHNICAL FRONT', "N/A", str(int(round(mech_f_xc,1))), str(int(round(mech_f_yc,1))), str(int(round(mech_f_R,1))), str(int(round(mech_f_d,1))), str(int(round(np.mean((mech_f_residu**2))**0.5,1)))],
	    ['MECHNICAL REAR', "N/A", str(int(round(mech_r_xc,1))), str(int(round(mech_r_yc,1))), str(int(round(mech_r_R,1))), str(int(round(mech_r_d,1))), str(int(round(np.mean((mech_r_residu**2))**0.5,1)))],
	    ['OPTICAL AXIS', str(round(np.abs(np.max([360*i['OA_reference_angle']/(2*np.pi) for i in OPT])), 2)), str(int(round(opt_xc,1))), str(int(round(opt_yc,1))), 
	    str(int(round(opt_R,1))), str(int(round(opt_d,1))), str(int(round(np.mean((opt_residu**2))**0.5,1)))]]
    print "[" + lens_configuration_name + "]"
    print 
    print tabulate(data, headers)
    print '\n'
    
  return MECH_F, MECH_R, OPT_F, OPT_R, OPT

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-pi", help="print info?", action='store_true')
  parser.add_argument("-p2d", help="plot 2d?", action='store_true')
  parser.add_argument("-p3d", help="plot 3d?", action='store_true')    
  args = parser.parse_args()
 
  CONFIG_FILE = "config.json"
  with open("config.json") as fp:
    cfg = json.load(fp)

  lens_configurations = np.array(cfg['PARAMS']['lens_configurations'])
  DB = cfg['PARAMS']['db_path']
  t_coordsys = cfg['TABLES']['T_COORDINATE_SYSTEMS']
  t_el = cfg['TABLES']['T_ELEMENTS']

  db = pymdb.parsefile(DB)
 
  for csys in cfg['COORDINATE_SYSTEMS']:
      '''
	  Process each coordinate system entry in the configuration file, 
	  cross-referencing it with the [t_coordsys] table on the [csMsRecNr] 
	  field and adding the appropriate transform as a new field.
      '''
      this_csys_csMsRecNr = csys['csMsRecNr'] 
      for entry in db[t_coordsys]['data']:	 
	  if this_csys_csMsRecNr == int(entry[db[t_coordsys]['headers'].index('csMsRecNr')]): # got it!
	      M11 = float(entry[db[t_coordsys]['headers'].index('csActM11')])
	      M12 = float(entry[db[t_coordsys]['headers'].index('csActM12')])
	      M13 = float(entry[db[t_coordsys]['headers'].index('csActM13')])
	      M14 = float(entry[db[t_coordsys]['headers'].index('csActM14')])
	      M21 = float(entry[db[t_coordsys]['headers'].index('csActM21')])
	      M22 = float(entry[db[t_coordsys]['headers'].index('csActM22')])
	      M23 = float(entry[db[t_coordsys]['headers'].index('csActM23')])
	      M24 = float(entry[db[t_coordsys]['headers'].index('csActM24')])
	      M31 = float(entry[db[t_coordsys]['headers'].index('csActM31')])
	      M32 = float(entry[db[t_coordsys]['headers'].index('csActM32')])
	      M33 = float(entry[db[t_coordsys]['headers'].index('csActM33')])
	      M34 = float(entry[db[t_coordsys]['headers'].index('csActM34')])
	      M41 = float(entry[db[t_coordsys]['headers'].index('csActM41')])
	      M42 = float(entry[db[t_coordsys]['headers'].index('csActM42')])
	      M43 = float(entry[db[t_coordsys]['headers'].index('csActM43')])
	      M44 = float(entry[db[t_coordsys]['headers'].index('csActM44')])

	      # actual values stored in database are the machine coordinates, so 
	      # the transform matrix is the coordinate transform required to convert to
	      # a specific coordinate system. for some reason, it's the inverse of the
	      # the transform matrix.
	      t = np.array([[M11, M12, M13, M14],
			  [M21, M22, M23, M24],
			  [M31, M32, M33, M34],
			  [M41, M42, M43, M44]])
	      t_inv = np.linalg.inv(t)

	      csys['t_inv'] = t_inv			# append it to the existing JSON object
	      
  for l in lens_configurations:
    solve(l, args)
    exit(0)
  

    

      





