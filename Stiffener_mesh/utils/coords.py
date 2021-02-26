import numpy as np

from utils.parameter import *
from utils.geometry import *

# ------- FLAT -------

def flat_point(param, dx, y):

	# dx range between 0 and param.dx (param.dx+1 points)
	x0  = 0
	discr  = param.dx

	length = 2.0 * param.X + param.fl

	x = (dx + 0.0) * length / discr + x0
	point = [x,y,0.0]

	return point

def limbflat_point_l(param, dx, y):

	# dx range between 0 and param.dl-1
	x0  = - param.limb
	length = param.limb
	discr  = param.dl

	x = (dx + 0.0) * length / discr + x0
	point = [x,y,0.0]

	return point

def limbflat_point_r(param, dx, y):

	# dx range between 1 and param.dl
	length = param.limb
	discr  = param.dl
	x0  = 2.0 * param.X + param.fl

	x = (dx + 0.0) * length / discr + x0
	point = [x,y,0.0]

	return point

def flat_point_side_r(param, dx, y):

	# dx range between 1 and param.dl
	length = param.side
	discr  = param.dl
	x0 = param.limb + 2.0 * param.X + param.fl

	x = (dx + 0.0) * length / discr + x0
	point = [x,y,0.0]

	return point

def flat_point_side_l(param, dx, y):

	# dx range between 0 and param.dl-1
	x0  = - param.side - param.limb
	length = param.side
	discr  = param.dl

	x = (dx + 0.0) * length / discr + x0
	point = [x,y,0.0]

	return point

# ------- CURVED -------

def limb_point_beg(param, dl, radius, Rotation):

	if param.Lshape==0:
		h = (-(param.limb) / (param.dl+0.0)) * dl + param.limb
		#h = (dl+0.0) * (param.limb) / (param.dl+0.0) * dl - param.limb
		#~ ([1*radius + 0*h,  0*radius + 1*h, 0.0]) # theta = 0
		#~ ([0*radius + 1*h, -1*radius + 0*h, 0.0]) # theta = np.pi/2
		point = np.asarray([radius*np.cos(Rotation) + h*np.sin(Rotation), -radius*np.sin(Rotation) + h*np.cos(Rotation), 0.0])
		# point = np.asarray([ h*np.sin(Rotation), h*np.cos(Rotation), 0.0])
	else:
		length = param.limb
		# length = param.limb + radius * (90.0-param.theta)*np.pi/180.0
		theta = param.theta * np.pi/180.
		h = length * (1- dl / param.dl)
		#~ ([ 0*radius + -1*h, -1*radius + 0*h, 0.0]) # theta = 0
		#~ ([-1*radius +  0*h,  0*radius + 1*h, 0.0]) # theta = np.pi/2
		point = np.asarray([h*np.cos(theta), h*np.sin(theta), 0.0])

	return point

def curve_point(param, dt, radius, Rotation):

	lame = 2.
	p = 2./lame

	theta = param.theta * np.pi/180.

	# t = np.power(dt/param.dc, 1./p) * (np.pi/2) + Rotation
	# point = np.asarray([radius*np.power(np.cos(t), p), -radius*np.power(np.sin(t), p), 0.0])

	if Rotation <0.001:
		x = np.sin(theta) * radius * np.power(1-dt/param.dc, 1./lame)
	else:
		x = - np.sin(theta) * radius * np.power(dt/param.dc, 1./lame)
	y = - radius * np.power(1 - np.power(np.abs(x), lame)/np.power(radius, lame), 1./lame)
	point = [x, y, 0.0]

	return point

def limb_point_end(param, dl, radius, Rotation):


	h = (dl+0.0)* (param.limb+0.0) / (param.dl+0.0) + 0.0
	#~ ([ 0*radius + -1*h, -1*radius + 0*h, 0.0]) # theta = 0
	#~ ([-1*radius +  0*h,  0*radius + 1*h, 0.0]) # theta = np.pi/2
	point = np.asarray([-radius*np.sin(Rotation) - h*np.cos(Rotation), -radius*np.cos(Rotation) + h*np.sin(Rotation), 0.0])

	return point
