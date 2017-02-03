import numpy as np

class optical_axis_solver():
  def __init__(self, lens_front_xyz, lens_rear_xyz):
    self.lens_front_xyz = np.array(lens_front_xyz)
    self.lens_rear_xyz = np.array(lens_rear_xyz)
    
  def _eval_direction_vector(self, normalise=True):
    dir_v = self.lens_front_xyz - self.lens_rear_xyz
    if normalise:
      return dir_v/np.linalg.norm(dir_v)
    else:
      return dir_v
    
  def getAngleBetweenOAAndDirectionVector(self, vector1, normalise_vector1=True, 
					  inDeg=False):
      '''
        Find angle from dot product = |dir_v_n|*|vector1_n|*cos(angle)
      '''
      dir_v_n = self._eval_direction_vector()
      vector1_n = vector1/np.linalg.norm(vector1)
      dotP = np.dot(dir_v_n, vector1)
      angle = np.arccos(dotP/(np.linalg.norm(dir_v_n)*np.linalg.norm(vector1))) 
      if inDeg:
	return 360*angle/(2*np.pi)
      else:
        return angle
          
  def getDirectionalCosines(self, inDeg=False):
    '''
      Find individual angles between optical axis direction vector, dir_v, and 
      component axes.
    '''
    dir_v_n = self._eval_direction_vector()
    cosX = dir_v_n[0]/np.sqrt((dir_v_n[0]**2)+(dir_v_n[1]**2)+(dir_v_n[2]**2))
    cosY = dir_v_n[1]/np.sqrt((dir_v_n[0]**2)+(dir_v_n[1]**2)+(dir_v_n[2]**2))
    cosZ = dir_v_n[2]/np.sqrt((dir_v_n[0]**2)+(dir_v_n[1]**2)+(dir_v_n[2]**2))
    angles = (np.arccos(cosX), np.arccos(cosY), np.arccos(cosZ))
    if inDeg:
      return 360*angles/2*np.pi
    else:
      return angles
  
  def getXY(self, z=0):
    '''
      Get (x, y) coordinates of optical axis at a given z.
    
      We're solving:
      
      x, y, z = lens_front_xyz(x, y, z) + t*dir_v(x, y, z)  ... (1)
      
      Reducing to parametric form for z,
      
      z = lens_front_xyz.z + t*dir_v.z ... (2)
    '''
    dir_v_n = self._eval_direction_vector()
    t = (z-self.lens_front_xyz[2])/dir_v_n[2]	# from (2)
    xyz_at_z = self.lens_front_xyz + t*dir_v_n	# from (1)
    
    return tuple(xyz_at_z[0:-1])