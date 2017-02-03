import pylab as plt
import numpy as np

class sag():
  def __init__(self, datasets, hard=False, fname="sag.png", inMicron=True):
    self.datasets = datasets
    self.hard = hard
    self.fname = fname
    self.inMicron = inMicron
 
  def _drawLinearDisplacementsToAxis(self, ax, x, y, x_err, y_err, angles, fit_xc, 
				     fit_yc, fit_r, title, color, lw1=3, lw2=1, 
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
      
    ax.errorbar(x, y, xerr=x_err, yerr=y_err, linewidth=lw1, label=title, color=color, 
		ls=ls1, elinewidth=elw1)
    for idx, angle in enumerate(angles):
      ax.annotate(int(round(360*angle/(2*np.pi))), xy=(x[idx], y[idx]), xytext=(-2, 1), 
		  textcoords='offset points', ha='right', va='bottom')
    ax.plot(fit_x, fit_y, linewidth=lw2, label='Circular fit', color=color, ls=ls2)
    ax.plot(fit_xc, fit_yc, color=color, marker=cmarker)
    
    ax.legend(loc='lower right', prop={'size':8})
    ax.set_title("Decentre relative to lens ring mount")
    ax.set_xlabel("x (" + unit + ")")
    ax.set_ylabel("y (" + unit + ")")
    ax.set(aspect=1, adjustable='datalim')
    
  def _drawRadialDisplacementsToAxis(self, ax, x, y, angles, title, color, lw=3, ls='-'):
    x = np.array(x)
    y = np.array(y)
    if self.inMicron:
      x *= 10**3
      y *= 10**3
      unit = 'micron'
    else:
      unit = 'mm'
      
    ax.plot(angles, [np.sqrt((x**2) + (y**2)) for x, y in zip(x, y)], linewidth=lw, 
	    color=color, ls=ls, label=title)
    plt.thetagrids(range(0,360,60), range(0,360,60))
    
  def _drawResidualsToAxis(self, ax, residuals, angles, title, color, lw=3, ls='-'):
    residuals = np.array(residuals)
    if self.inMicron:
      residuals *= 10**3
      unit = 'micron'
    else:
      unit = 'mm'
      
    ax.plot(angles, residuals, linewidth=lw, color=color, ls=ls, label=title)
    plt.thetagrids(range(0,360,60), range(0,360,60))
    ax.set_title("Residual from circular fit (" + unit + ")")
    
  def plot(self):
    '''
      This is a wrapper function to generate the complete sag plot.
      
      It requires datasets with two keys, 'data' and 'heading'. The former should 
      contain all necessary information (as a subdictionary) to call all the _draw* 
      functions.
    '''
    plot_colours = ('r', 'b', 'g', 'y')
    f, axes = plt.subplots(3, 1, figsize=(22,9))
    ax = plt.subplot(1, 3, 1)
    for idx, d in enumerate(self.datasets):
      self._drawLinearDisplacementsToAxis(ax, d['data']['x'], d['data']['y'], 
					  d['data']['x_err'], d['data']['y_err'], 
					  d['data']['angles'], d['data']['fit_xc'], 
					  d['data']['fit_yc'], d['data']['fit_r'],
					  d['heading'], 
					  color=plot_colours[idx])
    ax = plt.subplot(1, 3, 2, projection='polar')
    for idx, d in enumerate(self.datasets):
      self._drawRadialDisplacementsToAxis(ax, d['data']['x'], d['data']['y'], 
					  d['data']['angles'], title=d['heading'], 
					  color=plot_colours[idx])
    ax = plt.subplot(1, 3, 3, projection='polar')
    for idx, d in enumerate(self.datasets):
      self._drawResidualsToAxis(ax, d['data']['residuals'], d['data']['angles'], 
				title=d['heading'], color=plot_colours[idx])
    if self.hard:
      plt.savefig(self.fname)
    else:
      plt.show()
  