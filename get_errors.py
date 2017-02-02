'''
  TODO: blurb
'''

import pylab as plt
import numpy as np

from database import CMM_access

db = CMM_access("/home/barnsley/ELT-PCS/mech/CMM_Databases/SWIFT.QrtMeasDb")

X = []
Y = []
Z = []
DIM = []
for elMsRecNr in [31, 32, 33, 34, 35]:
  for elRecNr in db.getElementsFromelMsRecNr(elMsRecNr):
    el = db.getElementFromelRecNr(elRecNr)
    if el['elId'] == 'SPH_1':
      x, y, z = db.transElActPos1IntoPCS(elRecNr, 3)
      X.append(x)
      Y.append(y)
      Z.append(z)
      DIM.append(el['elActDim1'])
      
print np.mean(X), np.std(X)
print np.mean(Y), np.std(Y)
print np.mean(Z), np.std(Z)
print np.mean(DIM), np.std(DIM)
  