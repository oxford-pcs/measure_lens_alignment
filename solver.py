import numpy as np
from numpy.linalg import eig, inv

class axis():
  def __init__(self, pt1_xyz, pt2_xyz, pt1_radius, pt2_radius):
    self.pt1_xyz = np.array(pt1_xyz)
    self.pt2_xyz = np.array(pt2_xyz)
    self.pt1_radius = np.array(pt1_radius)
    self.pt2_radius = np.array(pt2_radius)
      
  def _eval_direction_cosines(self, vector):
    cosX = vector[0]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))
    cosY = vector[1]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))
    cosZ = vector[2]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))   
    
    return np.array([cosX, cosY, cosZ])
  
  def _getXYZAtLengthGivenOriginAndDirectionVector(self, origin, length, directionVector):
    '''
      Get (x, y, z) coordinates of point with length along axis.
      
      Consider,
      
      XYZ = XYZ(origin) + constant*directionVector ... (1)
      
      Moving to parametric form gives, e.g.
      
      X = X(origin) + constant*directionVector(x) ... (2)
      
      and knowing that,
      
      length =  sqrt(((X-X(origin))^2) + ((Y-Y(origin))^2) + ((Z-Z(origin))^2)) ... (3)
      
      We can substitute (3) into (2), solve for c and find XYZ by subsitution of c into (1)
    '''
    c = np.sqrt((length**2)/np.sum(directionVector**2))  
    return origin + (c*directionVector)
     
  def _eval_direction_vector(self, normalise=True, reverse=False):
    if np.isclose(self.pt1_xyz[2], self.pt2_xyz[2]):
      print "Element coordinates at same z."
      exit(0)
    
    if reverse:
      dir_v = self.pt2_xyz - self.pt1_xyz
    else:
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
    
  def getLensCentreThickness(self):
    '''
      Find lens centre thickness.
    '''
    if self.pt1_radius is None or self.pt2_radius is None:
      return "N/A"
    
    if self.pt1_xyz[2] > 0 and self.pt2_xyz[2] < 0: # biconcave
      XYZ_at_radius1 = self._getXYZAtLengthGivenOriginAndDirectionVector(self.pt1_xyz, self.pt1_radius, self._eval_direction_vector(reverse=True))
      XYZ_at_radius2 = self._getXYZAtLengthGivenOriginAndDirectionVector(self.pt2_xyz, self.pt2_radius, self._eval_direction_vector())
    elif self.pt1_xyz[2] < 0 and self.pt2_xyz[2] > 0: # biconvex
      XYZ_at_radius1 = self._getXYZAtLengthGivenOriginAndDirectionVector(self.pt1_xyz, self.pt1_radius, self._eval_direction_vector(reverse=True))
      XYZ_at_radius2 = self._getXYZAtLengthGivenOriginAndDirectionVector(self.pt2_xyz, self.pt2_radius, self._eval_direction_vector())  
    elif self.pt1_xyz[2] < 0 and self.pt2_xyz[2] < 0: # convex-concave
      XYZ_at_radius1 = self._getXYZAtLengthGivenOriginAndDirectionVector(self.pt1_xyz, self.pt1_radius, self._eval_direction_vector())
      XYZ_at_radius2 = self._getXYZAtLengthGivenOriginAndDirectionVector(self.pt2_xyz, self.pt2_radius, self._eval_direction_vector())        
    elif self.pt1_xyz[2] > 0 and self.pt2_xyz[2] > 0: # concave-convex
      XYZ_at_radius1 = self._getXYZAtLengthGivenOriginAndDirectionVector(self.pt1_xyz, self.pt1_radius, self._eval_direction_vector(reverse=True))
      XYZ_at_radius2 = self._getXYZAtLengthGivenOriginAndDirectionVector(self.pt2_xyz, self.pt2_radius, self._eval_direction_vector(reverse=True))  
    thickness = np.linalg.norm(XYZ_at_radius1 - XYZ_at_radius2)
    return thickness
    
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
  

        
