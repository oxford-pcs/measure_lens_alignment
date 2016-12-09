#!/usr/local/bin/python

# note z=0 is the origin of the reference surface for z, commonly the lens ring
import pymdb
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from tabulate import tabulate

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


CONFIG_FILE = "config.json"

with open("config.json") as fp:
  cfg = json.load(fp)

lenses = np.array(cfg['PARAMS']['lenses'])
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

nplots = 5
f, axes = plt.subplots(nplots, len(lenses))
for idx_l, l in enumerate(lenses):
    ''' 
        Now we process each lens configuration in turn.
    '''
    MECH_F = []		# front mechanical lens diameter
    MECH_R = []		# rear mechanical lens diameter
    OPT_F = []		# front lens surface
    OPT_R = []		# rear lens surface
    MECH = []		# front mechanical lens cylinder
    OPT = []		# optical axis
    this_lens_data = None
    for d in cfg['DATA']:					# which lens configuration are we processing?
        if d['id'] == l:
            this_lens_data = d
    if this_lens_data == None:      				# we couldn't find the configuration
        continue

      
    this_elMsRecNr_range = this_lens_data['elMsRecNr_range']	# get the range of elMsRecNr that correspond to this set of measurements
    for entry in db[t_el]['data']:
        this_elId = entry[db[t_el]['headers'].index('elId')]
        this_elMsRecNr = int(entry[db[t_el]['headers'].index('elMsRecNr')])
        for el in this_lens_data['elements']:
            if el['elId'] == this_elId and this_elMsRecNr in range(this_elMsRecNr_range[0], this_elMsRecNr_range[1]+1):
                for csys in cfg['COORDINATE_SYSTEMS']:		# find the corresponding coordinate system
                    if csys['csMsRecNr'] == this_elMsRecNr:
                        try:
                            t_inv = csys['t_inv']
                            x1 = float(entry[db[t_el]['headers'].index('elActPos1X')]) #elActPos1* are evaluated at z=0
                            y1 = float(entry[db[t_el]['headers'].index('elActPos1Y')])
                            z1 = float(entry[db[t_el]['headers'].index('elActPos1Z')])
                            x1y1z1_column_vector = np.array([[x1], [y1], [z1], [1]])
                            x1y1z1_transformed = np.dot(t_inv, x1y1z1_column_vector) 

                            elActDim1 = float(entry[db[t_el]['headers'].index('elActDim1')])
                            
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
    for f, r in zip(OPT_F, OPT_R):			
        if f['theta'] != r['theta']:		# should always be in same angle order, but check anyway
            exit(0)
        xyz_f = np.array([f['x'], f['y'], f['z']])
        xyz_r = np.array([r['x'], r['y'], r['z']])
        dir_v = xyz_f - xyz_r
        zis0 = -xyz_f[2]/dir_v[2]
        oa_at_zis0 = xyz_f + zis0*dir_v
        x = oa_at_zis0[0]
        y = oa_at_zis0[1]
        z = oa_at_zis0[2]
        r = np.sqrt((x**2) + (y**2) + (z**2))
        
        xyz_lens_surface_axis = np.array([0, 0, f['z']])
        
        dir_v_n = dir_v/np.linalg.norm(dir_v)
        xyz_lens_surface_axis_n = xyz_lens_surface_axis/np.linalg.norm(xyz_lens_surface_axis)
        
        OA_angle_lens_surface = np.arccos(np.dot(dir_v_n, xyz_lens_surface_axis_n))
        
        OPT.append({'x': x, 'y': y, 'z': z, 'theta': f['theta'], 'r': r, 'OA_reference_angle': OA_angle_lens_surface})
            
    '''
      Produce axis plots (mech and opt) in both linear and polar coordinates
    '''
    ax = plt.subplot(len(lenses), nplots, (idx_l*nplots)+1)
    
    # LINEAR
    ## MECHANICAL FRONT
    ### x, y
    x, y = [i['x'] for i in MECH_F], [i['y'] for i in MECH_F]
    ax.plot(x, y, linewidth=3, label='MECHANICAL FRONT', color='r', ls='-')
    for i in MECH_F:
        ax.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    ### circle
    mech_f_xc, mech_f_yc, mech_f_R, mech_f_residu = leastsq_circle([i['x'] for i in MECH_F], [i['y'] for i in MECH_F])
    theta_fit = np.linspace(-np.pi, np.pi, 180)
    mech_f_x_fit = mech_f_xc + mech_f_R*np.cos(theta_fit)
    mech_f_y_fit = mech_f_yc + mech_f_R*np.sin(theta_fit)
    ax.plot(mech_f_x_fit, mech_f_y_fit, color='r', ls='--', label='Circular fit')
    ax.plot(mech_f_xc, mech_f_yc, color='r', marker='x')
    
    mech_f_xd = x[-1]-x[0]
    mech_f_yd = y[-1]-y[0]
    mech_f_d = np.sqrt((mech_f_xd**2)+(mech_f_yd**2))

    ## MECHANICAL REAR
    ### x, y
    x, y = [i['x'] for i in MECH_R], [i['y'] for i in MECH_R]
    ax.plot(x, y, linewidth=3, label='MECHANICAL REAR', color='b', ls='-') 
    for i in MECH_R:
        ax.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    ### circle    
    mech_r_xc, mech_r_yc, mech_r_R, mech_r_residu = leastsq_circle([i['x'] for i in MECH_R], [i['y'] for i in MECH_R])
    theta_fit = np.linspace(-np.pi, np.pi, 180)
    mech_r_x_fit = mech_r_xc + mech_r_R*np.cos(theta_fit)
    mech_r_y_fit = mech_r_yc + mech_r_R*np.sin(theta_fit)
    ax.plot(mech_r_x_fit, mech_r_y_fit, color='b', ls='--', label='Circular fit')
    ax.plot(mech_r_xc, mech_r_yc, color='g', marker='x')
      
    mech_r_xd = x[-1]-x[0]
    mech_r_yd = y[-1]-y[0]
    mech_r_d = np.sqrt((mech_r_xd**2)+(mech_r_yd**2))
        
    ## OPTICAL AXIS  
    ### x, y
    x, y = [i['x'] for i in OPT], [i['y'] for i in OPT]
    ax.plot(x, y, linewidth=3, label='OPTICAL AXIS', color='g', ls='-')
    for i in OPT:
        ax.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    ### circle
    opt_xc, opt_yc, opt_R, opt_residu = leastsq_circle([i['x'] for i in OPT], [i['y'] for i in OPT])
    theta_fit = np.linspace(-np.pi, np.pi, 180)
    opt_x_fit = opt_xc + opt_R*np.cos(theta_fit)
    opt_y_fit = opt_yc + opt_R*np.sin(theta_fit)
    ax.plot(opt_x_fit, opt_y_fit, color='g', ls='--', label='Circular fit')
    ax.plot(opt_xc, opt_yc, color='b', marker='x')
  
    opt_xd = x[-1]-x[0]
    opt_yd = y[-1]-y[0]
    opt_d = np.sqrt((opt_xd**2)+(opt_yd**2))
    
    ax.legend(loc='lower right', prop={'size':8})
    ax.set_title("decentre relative to lens ring mount")
    ax.set_xlabel("x (micron)")
    ax.set_ylabel("y (micron)")
    ax.set(aspect=1, adjustable='datalim')

    # RADIAL
    ax = plt.subplot(len(lenses), nplots, (idx_l*nplots)+2, projection='polar')
    ax.plot([i['theta'] for i in MECH_F], [i['r'] for i in MECH_F], linewidth=3, color='r', ls='-', label='MECHANICAL FRONT')
    ax.plot([i['theta'] for i in MECH_R], [i['r'] for i in MECH_R], linewidth=3, color='g', ls='-', label='MECHANICAL REAR')
    ax.plot([i['theta'] for i in OPT], [i['r'] for i in OPT], linewidth=3, color='b', ls='-', label='OPTICAL AXIS')
    
    # RESIDUALS
    ax = plt.subplot(len(lenses), nplots, (idx_l*nplots)+3)
    ax.plot([int(round(360*i['theta']/(2*np.pi))) for i in MECH_F], mech_f_residu, color='r', ls='-', label='MECHANICAL FRONT')
    ax.plot([int(round(360*i['theta']/(2*np.pi))) for i in MECH_R], mech_r_residu, color='g', ls='-', label='MECHANICAL REAR')
    ax.plot([int(round(360*i['theta']/(2*np.pi))) for i in OPT], opt_residu, color='b', ls='-', label='OPTICAL AXIS')
    ax.set_xlim([0,360])
    ax.set_title("residual from circular fit")
    ax.set_xlabel("rotation angle (deg)")
    ax.set_ylabel("residual (micron)")
    
    headers = ['Axis', 'AXIS MAX DEVIATION (deg)', 'FIT_XCENT (micron)', 'FIT_YCENT (micron)', 'FIT_RADIUS (micron)', 'ELASTICITY (micron)', 'FIT_RESIDUAL_RMS (micron)']
    data = [['MECHNICAL FRONT', "N/A", str(int(round(mech_f_xc,1))), str(int(round(mech_f_yc,1))), str(int(round(mech_f_R,1))), str(int(round(mech_f_d,1))), str(int(round(np.mean((mech_f_residu**2))**0.5,1)))],
	    ['MECHNICAL REAR', "N/A", str(int(round(mech_r_xc,1))), str(int(round(mech_r_yc,1))), str(int(round(mech_r_R,1))), str(int(round(mech_r_d,1))), str(int(round(np.mean((mech_r_residu**2))**0.5,1)))],
	    ['OPTICAL AXIS', str(round(np.abs(np.max([360*i['OA_reference_angle']/(2*np.pi) for i in OPT])-np.min([360*i['OA_reference_angle']/(2*np.pi) for i in OPT])), 2)), str(int(round(opt_xc,1))), str(int(round(opt_yc,1))), 
            str(int(round(opt_R,1))), str(int(round(opt_d,1))), str(int(round(np.mean((opt_residu**2))**0.5,1)))]]
	     
    # RECESSION OF OPTICAL AXIS	 
    ax = plt.subplot(len(lenses), nplots, (idx_l*nplots)+4)
    ax.plot([int(round(360*i['theta']/(2*np.pi))) for i in OPT], [360*i['OA_reference_angle']/(2*np.pi) for i in OPT], 'x', color='b', label='OPTICAL AXIS')
    ax.set_xlim([0,360])
    ax.set_title("angle between OA and axis perpendicular\n to lens ring reference surface")
    ax.set_xlabel("rotation angle (deg)")
    ax.set_ylabel("angle (deg)")
    
    # RECESSION OF OPTICAL AXIS	 
    ax = plt.subplot(len(lenses), nplots, (idx_l*nplots)+5, projection='polar')
    ax.plot([i['theta'] for i in OPT], [360*i['OA_reference_angle']/(2*np.pi) for i in OPT], color='b', ls='-', label='OPTICAL AXIS')
    
    '''
      Pretty print pertinent data
    '''
    print "[" + l + "]"
    print 
    print tabulate(data, headers)
    print '\n'

plt.show()

    

     





