import numpy as np
from numpy.linalg import eig, inv
import transforms3d

class axis():
  def __init__(self, pt1_xyz, pt2_xyz, pt1_radius=None, pt2_radius=None, flip_lens=False, flip_PCS_z_direction=True, mount_ring_thickness=5):
    self.pt1_xyz = np.array(pt1_xyz)		# leftmost lens
    self.pt2_xyz = np.array(pt2_xyz)		# rightmost lens
    self.pt1_xyz[2]-=mount_ring_thickness	# move z origin so 0 is at the centre of the mount ring
    self.pt2_xyz[2]-=mount_ring_thickness
    if flip_lens:				# if lens has been measured in opposite orientation to how it is used
      self.pt1_xyz[2]*=-1
      self.pt2_xyz[2]*=-1  
      self.pt1_radius = np.array(pt2_radius)	
      self.pt2_radius = np.array(pt1_radius)
    else:
      self.pt1_radius = np.array(pt1_radius)	
      self.pt2_radius = np.array(pt2_radius)      

    if flip_PCS_z_direction:			# if the PCS z direction opposes the convention 
      self.pt1_xyz[2]*=-1
      self.pt2_xyz[2]*=-1
      
  def _eval_direction_cosines(self, vector):
    cosX = vector[0]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))
    cosY = vector[1]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))
    cosZ = vector[2]/np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))   
    
    return np.array([cosX, cosY, cosZ])
     
  def _eval_direction_vector(self, normalise=True, reverse=False):
    if np.isclose(self.pt1_xyz[2], self.pt2_xyz[2]):
      print "Element coordinates at same z."
      exit(0)

    # always head in positive Z
    if self.pt2_xyz[2]>self.pt1_xyz[2]:
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
      coordinate system axes.
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
    
    #
    # Since the direction vector is always evaluated in the positive Z direction, depending on the geometry of the lens, it may be 
    # necessary to flip this vector, e.g. think convex-concave
    #
    if self.pt1_xyz[2] > 0:
      XYZ_at_radius1 = self.getXYZAtLengthGivenOriginAndDirectionVector(self.pt1_xyz, self.pt1_radius, -1*self._eval_direction_vector())
    else:
      XYZ_at_radius1 = self.getXYZAtLengthGivenOriginAndDirectionVector(self.pt1_xyz, self.pt1_radius, self._eval_direction_vector())
    if self.pt2_xyz[2] > 0:
      XYZ_at_radius2 = self.getXYZAtLengthGivenOriginAndDirectionVector(self.pt2_xyz, self.pt2_radius, -1*self._eval_direction_vector())
    else:
      XYZ_at_radius2 = self.getXYZAtLengthGivenOriginAndDirectionVector(self.pt2_xyz, self.pt2_radius, self._eval_direction_vector())
    thickness = np.linalg.norm(XYZ_at_radius1 - XYZ_at_radius2)
    return thickness
  
  def getProjectedAngleInXYPlane(self, z=0, ref_axis=[0,1], centre=[0,0], inDeg=True):	
    '''
      Project the OA vector to z=z, calculate the XY position, construct a 
      2D vector from [centre] to this XY and measure the angle subtended by 
      this vector from [ref_axis] (clockwise).
    '''
    ref_axis = np.array(ref_axis)
    centre = np.array(centre)
    
    point_vector_from_fit_centre = np.array(self.getXY(z=z)) - centre
    dotP = np.dot(ref_axis, point_vector_from_fit_centre)
    crossP = np.cross(ref_axis, point_vector_from_fit_centre)
    angle = np.arccos(dotP/(np.linalg.norm(ref_axis)*np.linalg.norm(point_vector_from_fit_centre)))
    
    if np.sign(crossP) > 0:
      angle = (np.pi-angle) + np.pi
    
    if inDeg:
      return np.degrees(angle)
    else:
      return angle
  
  def getEulerAngles(self, align_axis=[0,0,1], order='xyz', rotating=True):
    '''
      Calculate the individual rotations required to align the OA vector with
      [align_axis]. Order is XYZ with an intrinsic rotating frame, this 
      is chosen to match that used by Zemax.
    '''
    align_axis = np.array(align_axis, dtype=float)
    dir_v = self._eval_direction_vector()

    dotP = np.dot(align_axis, dir_v)
    crossP = np.cross(align_axis, dir_v)
    angle = np.arccos(dotP)

    # this defines if the rotation axes are static or rotating
    if rotating:
      order = 'r' + order
    else:
      order = 's' + order

    return transforms3d.euler.axangle2euler(crossP, angle, order)
    
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
  
  def getXYZAtLengthGivenOriginAndDirectionVector(self, origin, length, directionVector):
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
  