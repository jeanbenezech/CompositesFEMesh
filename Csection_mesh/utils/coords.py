import numpy as np

from utils.parameters import *
from utils.geometry import *


def curve_point(param, dt, rep, stype):
	lame = 2.
	# p = 2./lame

	if param.substr[stype].Ori>0:
		theta =  np.pi/2.0
		radius = param.substr[stype].L-param.Y + param.substr[stype].thick
	else:
		theta =  -np.pi/2.0
		radius = param.substr[stype].L - param.substr[stype].thick

	center = np.asarray([radius*np.cos(theta), radius*np.sin(theta), 0.0])

	# t = np.power(dt/param.dc, 1./p) * (np.pi/2) + Rotation
	# point = np.asarray([radius*np.power(np.cos(t), p), -radius*np.power(np.sin(t), p), 0.0])

	if param.substr[stype].Ori>0:
		x = np.sin(param.substr[stype].Ori) * radius * np.power(dt/param.dc, 1./lame)
		y = - radius * np.power(1 - np.power(np.abs(x), lame)/np.power(radius, lame), 1./lame)
	else:
		x = - np.sin(param.substr[stype].Ori) * radius * np.power(dt/param.dc, 1./lame)
		y = radius * np.power(1 - np.power(np.abs(x), lame)/np.power(radius, lame), 1./lame)
	x+=center[0]
	y+=center[1]	
	
	point = [x*np.cos(rep.R)-y*np.sin(rep.R), y*np.cos(rep.R)+x*np.sin(rep.R), 0.0]
	new_center = [center[0]*np.cos(rep.R)-center[1]*np.sin(rep.R), center[1]*np.cos(rep.R)+center[0]*np.sin(rep.R), 0.0]
	# point = [x*np.cos(rep.R), y*np.sin(rep.R), 0.0]

	return point, new_center

def translation_along_line(substr):

	point = np.asarray([substr.L*np.cos(substr.Ori), substr.L*np.sin(substr.Ori), 0.0])
	return point

# def translation_along_radius(point, delta, Rotation):
# 	point = np.asarray([delta*np.sin(Rotation), -delta*np.cos(Rotation), 0.0])
