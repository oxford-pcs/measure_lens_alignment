import numpy as np
import datetime

import pymdb

class CMM_access():
  '''
      This class exposes routines to interact programatically with the 
      access database produced by a Wenzel CMM, and allows the user to 
      convert between the machine coordinate system (MCS) and part 
      coordinate systems (PCS).
      
      This program processes "measurements" in a Quartis CMM database. 
      Each measurement contains a group of "elements"; these are 
      circles, spheres, lines etc. 
      
      The Quartis CMM database has two pertinent tables to record the 
      information taken by the CMM: tbCoordSys and tbElement. The former 
      table contains the matrix elements of a transform that can be used 
      to convert the raw CMM table positions e.g. _tbElement.elAct* into 
      a desired PCS. The two tables can be joined on:
      
        _tbElement.elMsRecNr --> _tbCoordSys.csMsRecNr
           
      giving each element within a measurement a PCS (many elements will 
      thus have the same PCS).
      
      There may also exist several PCSs for each measurement. In this case,
      the user is required to specify a "csId", which uniquely identifies 
      a PCS in the database. The default for this is 1.
      
      For each element in a measurement, there will be a unique "elRecNr". 
  '''
  def __init__(self, db_file):
    self.db_file = db_file
    
    self.COORDSYS = []
    self.ELEMENTS = []
    self._parse()
    
  def _parse(self):
    '''
      Parse the database coordinate systems and elements into two separate 
      lists of dicts.
      
      Parsing the coordinate systems also adds the transform matrix which 
      would need to be applied to move from machine coordinates to the 
      corresponding PCS.
    '''
    db = pymdb.parsefile(self.db_file)
    
    cs_headers = db['_tbCoordSys']['headers']
    for cs_data in db['_tbCoordSys']['data']:
      csActM11 = float(cs_data[cs_headers.index('csActM11')])
      csActM12 = float(cs_data[cs_headers.index('csActM12')])
      csActM13 = float(cs_data[cs_headers.index('csActM13')])
      csActM14 = float(cs_data[cs_headers.index('csActM14')])
      csActM21 = float(cs_data[cs_headers.index('csActM21')])
      csActM22 = float(cs_data[cs_headers.index('csActM22')])
      csActM23 = float(cs_data[cs_headers.index('csActM23')])
      csActM24 = float(cs_data[cs_headers.index('csActM24')])
      csActM31 = float(cs_data[cs_headers.index('csActM31')])
      csActM32 = float(cs_data[cs_headers.index('csActM32')])
      csActM33 = float(cs_data[cs_headers.index('csActM33')])
      csActM34 = float(cs_data[cs_headers.index('csActM34')])
      csActM41 = float(cs_data[cs_headers.index('csActM41')])
      csActM42 = float(cs_data[cs_headers.index('csActM42')])
      csActM43 = float(cs_data[cs_headers.index('csActM43')])
      csActM44 = float(cs_data[cs_headers.index('csActM44')])
      
      # Actual values stored in database are the machine coordinates, so 
      # the 4x4 affine transformation matrix is the coordinate transform
      # required to convert to a specific coordinate system. For some  
      # reason, it's the inverse of the the transformation matrix.
      t = np.array([[csActM11, csActM12, csActM13, csActM14],
	            [csActM21, csActM22, csActM23, csActM24],
		    [csActM31, csActM32, csActM33, csActM34],
		    [csActM41, csActM42, csActM43, csActM44]])
      t_inv = np.linalg.inv(t)
      	    
      self.COORDSYS.append({'csMsRecNr': int(cs_data[cs_headers.index('csMsRecNr')]),
			    'csId': int(cs_data[cs_headers.index('csId')]),			# may have multiple PCS
			    'csActM11' : csActM11,
			    'csActM12' : csActM12,
			    'csActM13' : csActM13,
			    'csActM14' : csActM14,
			    'csActM21' : csActM21,
			    'csActM22' : csActM22,
			    'csActM23' : csActM23,
			    'csActM24' : csActM24,
			    'csActM31' : csActM31,
			    'csActM32' : csActM32,
			    'csActM33' : csActM33,
			    'csActM34' : csActM34,
			    'csActM41' : csActM41,
			    'csActM42' : csActM42,
			    'csActM43' : csActM43,
			    'csActM44' : csActM44,  
			    't': t,
	                    't_inv': t_inv
                           })
      
    el_headers = db['_tbElement']['headers']
    for el_data in  db['_tbElement']['data']:
      self.ELEMENTS.append({'elRecNr' : int(el_data[el_headers.index('elRecNr')]),
			    'elMsRecNr' : int(el_data[el_headers.index('elMsRecNr')]),
			    'elId' : str(el_data[el_headers.index('elId')]),			# CIR_1, PLN_1 etc.
			    'elTime': datetime.datetime.strptime(str(el_data[el_headers.index('elTime')]), '%Y.%m.%d %H:%M:%S:%f'),
			    'elType' : int(el_data[el_headers.index('elType')]),
			    'elActPos1X' : float(el_data[el_headers.index('elActPos1X')]),	# machine coordinate, projected to a reference surface (x)
			    'elActPos1Y' : float(el_data[el_headers.index('elActPos1Y')]),	# machine coordinate, projected to a reference surface (y)
			    'elActPos1Z' : float(el_data[el_headers.index('elActPos1Z')]),	# machine coordinate, projected to a reference surface (z)
			    'elActPos2X' : float(el_data[el_headers.index('elActPos2X')]),      # machine coordinate
			    'elActPos2Y' : float(el_data[el_headers.index('elActPos2Y')]),      # machine coordinate
			    'elActPos2Z' : float(el_data[el_headers.index('elActPos2Z')]),      # machine coordinate
			    'elActDir1X' : float(el_data[el_headers.index('elActDir1X')]),	# direction vector (x)
			    'elActDir1Y' : float(el_data[el_headers.index('elActDir1Y')]),	# direction vector (y)
			    'elActDir1Z' : float(el_data[el_headers.index('elActDir1Z')]),	# direction vector (z)
			    'elActDir2X' : float(el_data[el_headers.index('elActDir2X')]),	# blank (not projected!)
			    'elActDir2Y' : float(el_data[el_headers.index('elActDir2Y')]),      # blank (not projected!)
			    'elActDir2Z' : float(el_data[el_headers.index('elActDir2Z')]),      # blank (not projected!)
			    'elActDim1' : float(el_data[el_headers.index('elActDim1')]),
			    'elActDim2' : float(el_data[el_headers.index('elActDim2')]),	
			    'elRefKind' : str(el_data[el_headers.index('elRefKind')])		# reference element
			    })
 
  def getElementsFromelMsRecNr(self, elMsRecNr):
    '''
      Get elRecNr for all elements with a given elMsRecNr. 
      
      If an elId is specified, the results will be filtered to show only that entry
      with the given elId.
    '''
    entries = [entry['elRecNr'] for entry in self.ELEMENTS if entry['elMsRecNr'] == elMsRecNr]
    try:
      assert len(entries) > 0
      return entries
    except AssertionError:
      print "No entries found with given elMsRecNr"
      
  def getElementFromelRecNr(self, elRecNr):
    '''
      Get element from its elRecNr.
    '''
    try:
      el_entry = (entry for entry in self.ELEMENTS if entry['elRecNr'] == elRecNr).next()
      return el_entry 
     
    except StopIteration:
      print "elRecNr not found."   
      
  def transElActPos1IntoPCS(self, elRecNr, csId=1, projected=False):
    '''
      Transform an element's elActPos1* values (evauluated at z=0) into the PCS.
    '''
    try:
      el_entry = self.getElementFromelRecNr(elRecNr)
      cs_entry = (entry for entry in self.COORDSYS if entry['csMsRecNr'] == el_entry['elMsRecNr'] and entry['csId'] == csId).next()
      
      if not projected:
        x1 = el_entry['elActPos2X']
        y1 = el_entry['elActPos2Y']
        z1 = el_entry['elActPos2Z']
      else:
        x1 = el_entry['elActPos1X']
        y1 = el_entry['elActPos1Y']
        z1 = el_entry['elActPos1Z']
      x1y1z1_column_vector = np.array([[x1], [y1], [z1], [1]])	
      x1y1z1_transformed = np.dot(cs_entry['t_inv'], x1y1z1_column_vector)
      
      return x1y1z1_transformed[:-1].flatten()
    
    except StopIteration:
      print "Couldn't find coordinate system."
     
      

      
