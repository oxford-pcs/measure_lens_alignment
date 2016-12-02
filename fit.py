import numpy as np
import pylab as plt
from scipy      import optimize

def fit_sphere(x, y, z):
  # see http://jekel.me/2015/Least-Squares-Sphere-Fit/
  #   Assemble the A matrix
  x = np.array(x)
  y = np.array(y)
  z = np.array(z)
  A = np.zeros((len(x),4))
  A[:,0] = x*2
  A[:,1] = y*2
  A[:,2] = z*2
  A[:,3] = 1
  
  #   Assemble the f matrix
  f = np.zeros((len(x),1))
  f[:,0] = (x*x) + (y*y) + (z*z)
  C, residuals, rank, singval = np.linalg.lstsq(A,f)

  #   solve for the radius
  t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
  radius = np.sqrt(t)
  
  return radius[0], np.array([C[0][0], C[1][0], C[2][0]])

def calculate_3d_pt_pt_vector(v1, v2):
  '''
    x, y, z = v1(x, y, z) + t*v_dir(x, y, z)
    
    parametric form:
    
    z = v1_z + t*v_dir_z
    solve for z=0:
    
    t = -v1_z/v_dir_z
  '''
  v_dir = v2 - v1
  zis0 = -v1[2]/v_dir[2]
  return v1 + zis0*v_dir
  
def parse_file(fname):
  x = []
  y = []
  z = []
  with open(fname) as f:
    for line in f:
      try:
	x.append(float(line.split()[1]))
	y.append(float(line.split()[2]))
	z.append(float(line.split()[3]))
      except ValueError:
	pass
  return np.array(x), np.array(y), np.array(z)

if __name__ == "__main__":
  data = []
  with open("order") as f:
    for line in f:
      ftype = line.split()[0]
      if "opt" in ftype:
	f1 = line.split()[1]
	f2 = line.split()[2]
	x1, y1, z1 = parse_file(f1)
	x2, y2, z2 = parse_file(f2)	
	rad1, c1 = fit_sphere(x1, y1, z1)
	rad2, c2 = fit_sphere(x2, y2, z2)

	r = calculate_3d_pt_pt_vector(c1, c2)*10**3
	
	data.append({"COORDS": (r[0], r[1], r[2]), 
		    "TITLE": ftype,
		    "STYLE": 'bo'
		    })
	
data.append({"COORDS": (-0.005*10**3, 0.452*10**3, 0), 
	     "TITLE": "mech_1f",
	     "STYLE": "ro"
	     })
     	    
fig = plt.figure()
for entry in data:
  plt.plot(entry['COORDS'][0], entry['COORDS'][1], entry['STYLE'], label=entry['TITLE'])
plt.xlabel("displacement x (micron)")
plt.ylabel("displacement y (micron)")
plt.legend(numpoints=1)
plt.show()
    
