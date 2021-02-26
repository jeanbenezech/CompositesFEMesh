import numpy as np

from utils.parameter import *

def w_curve(param, s, r, z, sdef=0, rdef=0, zdef=0):
	
	if (sdef==0 and rdef==0 and zdef==0):
		sdef = param.limb + (param.X - param.Y/2.0)*np.pi/4.0
		#~ sdef = (param.X - param.Y/2.0)*np.pi/4.0
		rdef = param.X - param.Y/2.0
		#~ rdef = param.X - param.Y
		zdef = param.Z/2.0
	
	d=1.0
	b1=1.0/4.0
	b2=1.0/1.0
	b3=1.0/2.0
	S=2*param.limb+(param.X - param.Y/2.0)*np.pi/2.0
	R=param.X-param.Y # min radius
	W=param.Z
	
	#~ delta = d * 1./(np.cosh((s-sdef)/b1)**2) * \
	delta = d * 1./(np.cosh(np.pi*((s-sdef)/(S*b1)))**2) * \
				1./(np.cosh(np.pi*((r-rdef)/(R*b2)))**2) * \
				1./(np.cosh(np.pi*((z-zdef)/(W*b3)))**2)
	
	#~ print ('sdef, rdef, zdef : ', sdef, rdef, zdef)
	#~ print ('s, r, z : ', s, r, z)
	#~ print ('delta : ', delta)
	#~ print ()
	
	return delta

def w_flat(param, s, r, z, sdef=0, rdef=0, zdef=0):
	
	if (sdef==0 and rdef==0 and zdef==0):
		sdef = param.X/2.0
		rdef = -param.Y/2.0
		zdef = param.Z/2.0
		
	
	d=1.0
	b1=1.0/4.0
	b2=1.0/1.0
	b3=1.0/2.0
	S=param.X
	R=param.Y
	W=param.Z
	
	delta = d * 1./(np.cosh(np.pi*((s-sdef)/(S*b1)))**2) * \
				1./(np.cosh(np.pi*((r-rdef)/(R*b2)))**2) * \
				1./(np.cosh(np.pi*((z-zdef)/(W*b3)))**2)
	
	#~ print ('sdef, rdef, zdef : ', sdef, rdef, zdef)
	#~ print ('s, r, z : ', s, r, z)
	#~ print ('delta : ', delta)
	#~ print ()
	
	return delta
