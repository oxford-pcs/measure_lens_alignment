import pylab as plt
import numpy as np

class sag():
  def __init__(self, datasets, inMicron=True):
    self.datasets = datasets
    self.inMicron = inMicron
    
  def _drawAnglesFromMountNormalToAxis(self, ax, xy_angles_from_12_o_clock, angles_from_mount_normal, mount_angles, label, color, lw=3, ls='-'):
    xy_angles_from_12_o_clock = np.array(xy_angles_from_12_o_clock)
    angles_from_mount_normal = np.array(angles_from_mount_normal)*60

    for idx, angle in enumerate(mount_angles):
      ax.annotate(int(round(360*angle/(2*np.pi))), xy=(xy_angles_from_12_o_clock[idx], angles_from_mount_normal[idx]), xytext=(-2, 1), 
		  textcoords='offset points', ha='right', va='bottom', color=color, fontsize=14)
    ax.plot(xy_angles_from_12_o_clock, angles_from_mount_normal, linewidth=lw, color=color, ls=ls, label=label)
    ax.set_theta_offset(np.pi)
    ax.set_theta_direction(-1)
    plt.thetagrids(range(0,360,60), range(0,360,60))
    ax.set_title("Angle between axis and mount frame\nnormal [0, 0, 1] (arcmin)\n")
 
  def _drawLinearDisplacementsToAxis(self, ax, x, y, x_err, y_err, mount_angles, fit_xc, 
				     fit_yc, fit_r, label, color, lw1=3, lw2=1, 
				     ls1='-', ls2='--', cmarker='x', elw1=1):
    x = np.array(x)
    y = np.array(y)
    if self.inMicron:
      x *= 10**3
      y *= 10**3
      x_err *= 10**3
      y_err *= 10**3
      fit_xc *= 10**3
      fit_yc *= 10**3
      fit_r *= 10**3
      unit = 'micron'
    else:
      unit = 'mm'
      
    theta_fit = np.linspace(-np.pi, np.pi, 180)
    fit_x = fit_xc + fit_r*np.cos(theta_fit)
    fit_y = fit_yc + fit_r*np.sin(theta_fit)
      
    ax.errorbar(x, y, xerr=x_err, yerr=y_err, linewidth=lw1, label=label, color=color, 
		ls=ls1, elinewidth=elw1)
    for idx, angle in enumerate(mount_angles):
      ax.annotate(int(round(360*angle/(2*np.pi))), xy=(x[idx], y[idx]), xytext=(-2, 1), 
		  textcoords='offset points', ha='right', va='bottom', color=color, fontsize=14)
    ax.plot(fit_x, fit_y, linewidth=lw2, label='Circular fit', color=color, ls=ls2)
    ax.plot(fit_xc, fit_yc, color=color, marker=cmarker)
    
    ax.legend(loc='lower right', prop={'size':10})
    ax.set_title("xy (z=-mount_ring_thickness/2) decentre\nrelative to lens ring mount\n")
    ax.set_xlabel("x (" + unit + ")")
    ax.set_ylabel("y (" + unit + ")")
    ax.set(aspect=1, adjustable='datalim')
    
  def _drawRadialDisplacementsToAxis(self, ax, xy_angles_from_12_o_clock, xy, mount_angles, label, color, lw=3, ls='-'):
    xy = np.array(xy)
    if self.inMicron:
      xy *= 10**3
      unit = 'micron'
    else:
      unit = 'mm'
    x = xy[0]
    y = xy[1]
    radial_displacements = np.sqrt((x**2) + (y**2))
    
    for idx, angle in enumerate(mount_angles):
      ax.annotate(int(round(360*angle/(2*np.pi))), xy=(xy_angles_from_12_o_clock[idx], radial_displacements[idx]), 
		  xytext=(-2, 1), textcoords='offset points', ha='right', va='bottom', color=color, fontsize=14)

    ax.plot(xy_angles_from_12_o_clock, radial_displacements, linewidth=lw, 
	    color=color, ls=ls, label=label)
    plt.thetagrids(range(0,360,60), range(0,360,60))
    ax.set_theta_offset(np.pi)
    ax.set_theta_direction(-1)
    ax.set_title("Radial decentre from circular fit\n as a function of the angle between the fit centre\n and vector [0, 1], clockwise (" + unit + ")\n")
    
  def _drawResidualsToAxis(self, ax, xy_angles_from_12_o_clock, residuals, mount_angles, label, color, 
			   lw=3, ls='-'):
    residuals = np.array(residuals)
    if self.inMicron:
      residuals *= 10**3
      unit = 'micron'
    else:
      unit = 'mm'
      
    for idx, angle in enumerate(mount_angles):
      ax.annotate(int(round(360*angle/(2*np.pi))), xy=(xy_angles_from_12_o_clock[idx], residuals[idx]), 
		  xytext=(-2, 1), textcoords='offset points', ha='right', va='bottom', color=color, fontsize=14)  
      
    ax.plot(xy_angles_from_12_o_clock, residuals, linewidth=lw, color=color, ls=ls, label=label)
    plt.thetagrids(range(0,360,60), range(0,360,60))
    ax.set_theta_offset(np.pi)
    ax.set_theta_direction(-1)
    ax.set_title("Residual from circular fit (" + unit + ")\n")
    
  def draw(self, hard=False, fname="sag.png"):
    if hard:
      plt.savefig(fname)
    else:
      plt.show() 
    
  def plot(self):
    '''
      This is a wrapper function to generate the complete sag plot.
      
      It requires datasets with two keys, 'data' and 'heading'. The former should 
      contain all necessary information (as a subdictionary) to call all the _draw* 
      functions.
    '''
    plot_colours = ('r', 'b', 'g', 'y')
    f, axes = plt.subplots(3, 1, figsize=(16,7))
    ax = plt.subplot(1, 4, 1)
    plt.tick_params(labelsize=10)
    plt.rcParams.update({'axes.titlesize': 'small', 'axes.labelsize': 'small', 'xtick.labelsize':'small', 'ytick.labelsize':'small'})
    for idx, d in enumerate(self.datasets):
      self._drawLinearDisplacementsToAxis(ax, d['data']['x'], d['data']['y'], 
					  d['data']['x_err'], d['data']['y_err'], 
					  d['data']['mount_angles'], d['data']['fit_xc'], 
					  d['data']['fit_yc'], d['data']['fit_r'],
					  d['heading'], 
					  color=plot_colours[idx])
    ax = plt.subplot(1, 4, 2, projection='polar')
    for idx, d in enumerate(self.datasets):
      self._drawRadialDisplacementsToAxis(ax, d['data']['xy_angles_from_12_o_clock'],
					  (d['data']['x'], d['data']['y']), 
					  d['data']['mount_angles'], label=d['heading'], 
					  color=plot_colours[idx])
    ax = plt.subplot(1, 4, 3, projection='polar')
    for idx, d in enumerate(self.datasets):
      self._drawResidualsToAxis(ax, d['data']['xy_angles_from_12_o_clock'],
				d['data']['residuals'], d['data']['mount_angles'], 
				label=d['heading'], color=plot_colours[idx])
    ax = plt.subplot(1, 4, 4, projection='polar')
    for idx, d in enumerate(self.datasets):
      self._drawAnglesFromMountNormalToAxis(ax, d['data']['xy_angles_from_12_o_clock'],
					    [angle[2] for angle in 
					     d['data']['angles_from_mount_normal']],
                                            d['data']['mount_angles'],
                                            label=d['heading'], color=plot_colours[idx])
