import numpy as np
from numpy.linalg import eig, inv

class axis():
  def __init__(self, pt1_xyz, pt2_xyz):
    self.pt1_xyz = np.array(pt1_xyz)
    self.pt2_xyz = np.array(pt2_xyz)
      
  def _eval_direction_cosines(self, vector):
    cosX = vector[0]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))
    cosY = vector[1]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))
    cosZ = vector[2]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))   
    
    return np.array([cosX, cosY, cosZ])
    
  def _eval_direction_vector(self, normalise=True):
    if np.isclose(self.pt1_xyz[2], self.pt2_xyz[2]):
      print "Element coordinates at same z."
      exit(0)
    
    dir_v = self.pt1_xyz - self.pt2_xyz
    if normalise:
      return dir_v/np.linalg.norm(dir_v)
    else:
      return dir_v
    
  def getAngleBetweenOAAndDirectionVector(self, vector1, inDeg=False):
      '''
        Find angle from dot product = |dir_v|*|vector1|*cos(angle)
      '''
      dir_v = self._eval_direction_vector(normalise=False)
      dotP = np.dot(dir_v, vector1)
      angle = np.arccos(dotP/(np.linalg.norm(dir_v)*np.linalg.norm(vector1)))
       
      if dir_v[2] < 0: # direction vector z points in opposite direction to reference vector		
	angle = np.pi-angle
      if inDeg:
	return (360*angle)/(2*np.pi)
      else:
        return angle
          
  def getComponentAngles(self, inDeg=False):
    '''
      Find individual angles between optical axis direction vector, dir_v, and 
      component axes.
    '''
    dir_v_n = self._eval_direction_vector()
    angles = np.array(np.arccos(self._eval_direction_cosines(dir_v_n)))
    
    if dir_v_n[2] < 0: # direction vector z points in opposite direction to reference vector
      angles[2] = np.pi-angles[2]
    if inDeg:
      return (360*angles)/(2*np.pi)
    else:
      return angles
  
  def getXY(self, z=0):
    '''
      Get (x, y) coordinates of optical axis at a given z.
    
      We're solving:
      
      x, y, z = pt1_xyz(x, y, z) + t*dir_v(x, y, z)  ... (1)
      
      Reducing to parametric form for z,
      
      z = pt1_xyz.z + t*dir_v.z ... (2)
    '''
    dir_v_n = self._eval_direction_vector()
    t = (z-self.pt1_xyz[2])/dir_v_n[2]	# from (2)
    xyz_at_z = self.pt1_xyz + t*dir_v_n	# from (1)
    
    return tuple(xyz_at_z[0:-1])
  
        
