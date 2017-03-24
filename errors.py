import numpy as np
from scipy.spatial import distance

class measurementError():
  '''
    Requires repeated measurements at single position (no rotation).
    
    Takes a set of solver.axis instances as input.
  '''
  def __init__(self, axes):
    self.axes = axes
    
  def calculate_angle_error_at_z(self, z=0, vector=[0,0,1]):
    '''
      Returns the standard deviation of the set of angles between the OA and [vector].
    '''
    angular_deviation = []
    for ax in self.axes:
      angular_deviation.append(ax.getAngleBetweenOAAndDirectionVector(vector, inDeg=True))
    
    return np.std(angular_deviation)    
    
  def calculate_lens_thickness_error(self):
    return np.std([ax.getLensCentreThickness() for ax in self.axes])
  
  def calculate_lens_radius1_error(self):
    return np.std([ax.pt1_radius for ax in self.axes])  
  
  def calculate_lens_radius2_error(self):
    return np.std([ax.pt2_radius for ax in self.axes])    
  
  def calculate_position_error_at_z(self, z=0):
    '''
      Returns the standard deviation in x and y, and the euclidean distance between 
      pairs of coordinates.    
    '''
    xy_at_given_z = []
    for ax in self.axes:
      x, y = ax.getXY(z=z)
      xy_at_given_z.append((x,y))
    X = [xy[0] for xy in xy_at_given_z]
    Y = [xy[1] for xy in xy_at_given_z]
    
    pairs = []
    for x, y in zip(X, Y):
      pairs.append((x,y))
      
    distances = distance.pdist(pairs)

    return ((np.std(X), np.std(Y)), np.mean(distances))
  