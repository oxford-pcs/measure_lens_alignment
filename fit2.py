# note z=0 is the origin of the reference surface for z, commonly the lens ring
import pymdb
import json
import numpy as np
import matplotlib.pyplot as plt

CONFIG_FILE = "config.json"
DB_FILE = "SWIFT.mdb"

with open("config.json") as fp:
  cfg = json.load(fp)

lenses = np.array(cfg['PARAMS']['lenses'])

t_coordsys = cfg['TABLES']['T_COORDINATE_SYSTEMS']
t_el = cfg['TABLES']['T_ELEMENTS']

db = pymdb.parsefile(DB_FILE)

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

MECH_F = []
MECH_R = []
OPT_F = []
OPT_R = []
for l in lenses:
    ''' 
        Now we process each lens configuration in turn
    '''
    for d in cfg['DATA']:			# which lens configuration are we processing?
        this_lens_data = None
        if d['id'] == l:
            this_lens_data = d
    if this_lens_data == None:      # we couldn't find the configuration
        continue
    for entry in db[t_el]['data']:
        this_elId = entry[db[t_el]['headers'].index('elId')]
        this_elMsRecNr = int(entry[db[t_el]['headers'].index('elMsRecNr')])
        for el in this_lens_data['elements']:
            if el['elId'] == this_elId:
                for csys in cfg['COORDINATE_SYSTEMS']:	# find the corresponding coordinate system
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
    OPT = []
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
        OPT.append({'x': x, 'y': y, 'z': z, 'theta': f['theta'], 'r': r, 'radius': None})
            
    '''
      Produce axis plots (mech and opt) in both linear and polar coordinatesss
    '''
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax = plt.subplot(121)
    ax.set_xlim([-50, 50])
    ax.set_ylim([-50, 50])
    ax.plot([i['x'] for i in MECH_F], [i['y'] for i in MECH_F], 'o-', linewidth=3, label='MECHANICAL FRONT')
    for i in MECH_F:
        plt.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    ax.plot([i['x'] for i in MECH_R], [i['y'] for i in MECH_R], 'o-', linewidth=3, label='MECHANICAL REAR')
    for i in MECH_R:
        plt.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    ax.plot([i['x'] for i in OPT], [i['y'] for i in OPT], 'o-', linewidth=3, label='OPTICAL AXIS')
    for i in OPT:
        plt.annotate(int(round(360*i['theta']/(2*np.pi))), xy = (i['x'], i['y']), xytext = (-2, 1), textcoords = 'offset points', ha = 'right', va = 'bottom')
    plt.legend(loc='lower right')
    ax = plt.subplot(122, projection='polar')
    ax.plot([i['theta'] for i in MECH_F], [i['r'] for i in MECH_F], linewidth=3)
    ax.plot([i['theta'] for i in MECH_R], [i['r'] for i in MECH_R], linewidth=3)
    ax.plot([i['theta'] for i in OPT], [i['r'] for i in OPT], linewidth=3)
    plt.show()

    

     





