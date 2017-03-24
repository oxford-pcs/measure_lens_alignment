'''
    The prerequisite to these classes is that the user has obtained a list of 
    axis (either mechanical or optical) x, y coordinates (probably using 
    solver.py) for various lens rotations. 
'''

import numpy as np
from scipy import optimize

class sag():
  ''' 
    Requires measurements to have been made at each rotation angle.
    
    To calculate sag, a least squares circular fit is made to the list of 
    x, y coordinates, the residual of which gives an indication of the lens 
    sag.
    
    If the lens mount was used to define the PCS, the radius of the fit 
    will be equivalent to the displacement of the axis centre from the 
    mount centre.
    
    Takes a set of solver.axis instances as input.
  '''
  def __init__(self, axes):
    self.axes = axes

  def _optLeastSqCircle(self, x,y):
    def calcR(x,y, xc, yc):
      '''
        Calculate distance of each point from the center (xc, yc) .
      '''
      return np.sqrt((x-xc)**2 + (y-yc)**2)

    def f(c, x, y):
      '''
        Calculate the algebraic distance between the data points and the mean 
        circle centered at c = (xc, yc).
      '''
      Ri = calcR(x, y, *c)
      return Ri - Ri.mean()
    
    x_m = np.mean(x)
    y_m = np.mean(y)
    centre_estimate = x_m, y_m
    centre, ier = optimize.leastsq(f, centre_estimate, args=(x,y))
    xc, yc = centre
    Ri = calcR(x, y, *centre)
    R = Ri.mean()
    residuals = np.sqrt((Ri - R)**2)
    
    return xc, yc, R, residuals
	    
  def calculate(self, z=0):
    xy = []
    for ax in self.axes:
      xy.append(ax.getXY(z=z))
    x = [entry[0] for entry in xy]
    y = [entry[1] for entry in xy]
    xc, yc, R, residuals = self._optLeastSqCircle(x, y)
    return (xc, yc, R, residuals)
  
class hysteresis():
  '''
    Requires a repeated measurement of the initial angle, taken after a full 
    rotation has been made.
    
    To evaluate hysteresis, defined here as the ability of the lens to return 
    to its nominal position, the user can give two indices in the x, y list
    which (should) correspond to the same rotation angle. The difference in x, y 
    between these coordinates is then calculated.
  '''
  def __init__(self, axes, angles):
    self.axes = axes
    self.angles = angles
    
  def calculate(self, idx1, idx2, z):
    xy = []
    for ax in self.axes:
      xy.append(ax.getXY(z=z))
    x = [entry[0] for entry in xy]
    y = [entry[1] for entry in xy]
    try:
      xd = x[idx1]-x[idx2]
      yd = y[idx1]-y[idx2]
      d = np.sqrt((xd**2)+((yd**2))) 
    except IndexError:
      print "One of the indices is out of range."
      return None
    
    return d
  
  
