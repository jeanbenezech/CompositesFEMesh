import numpy as np

from utils.parameter import *
from utils.geometry import *
from utils.coords import *
from utils.gmsh_hierarchy import *

# ~~~~~~~~~~~~~~ ADD POINTS ~~~~~~~~~~~~~~

def Flat_sides(geo, param, cntP, cntL, Translation, y, dy=0):
	z = Translation[2]
	r = y

	# ~~~~~~~~ RIGHT ~~~~~~~~

	points=[]
	line=[]

	line.append(cntP-param.dl-1) # param.dl:: discr of limb
	for dx in range(1, param.dl+1):
		point = flat_point_side_r(param, dx, y)
		points.append(point + Translation)
		line.append(cntP)
		cntP+=1

	geo.All_points.append(points)
	cntL=add_x_line(geo, param, line, cntL, dy, bool_side=1)

	# ~~~~~~~~ LEFT ~~~~~~~~

	points=[]
	line=[]

	for dx in range(param.dl):
		point = flat_point_side_l(param, dx, y)
		points.append(point + Translation)
		line.append(cntP)
		cntP+=1

	line.append(cntP - (3*param.dl))

	geo.All_points.append(points)
	cntL=add_x_line(geo, param, line, cntL, dy, bool_side=1)

	return cntP, cntL

def Flat_points(geo, param, cntP, cntL, Translation, y, dy=0):

	z = Translation[2]
	r = y

	# ~~~~~~~~ MIDDLE ~~~~~~~~
	points=[]
	line=[]

	for dx in range(param.dx+1):
		point = flat_point(param, dx, y)
		points.append(point + Translation)
		line.append(cntP)
		cntP+=1

	geo.All_points.append(points)
	cntL=add_x_line(geo, param, line, cntL, dy, 1, 0)

	# ~~~~~~~~ RIGHT ~~~~~~~~
	points=[]
	line=[]

	line.append(cntP-1)
	for dl in range(1, param.dl+1):
		point = limbflat_point_r(param, dl, y)
		points.append(point + Translation)
		line.append(cntP)
		cntP+=1

	geo.All_points.append(points)
	cntL=add_x_line(geo, param, line, cntL, dy, 1)

	# ~~~~~~~~ LEFT ~~~~~~~~
	points=[]
	line=[]

	for dl in range(param.dl):
		point = limbflat_point_l(param, dl, y)
		points.append(point + Translation)
		line.append(cntP)
		cntP+=1

	line.append(cntP-(param.dx+1+2*param.dl))

	geo.All_points.append(points)
	cntL=add_x_line(geo, param, line, cntL, dy, 1)

	return cntP, cntL

def Lshape_points(geo, param, radius, cntP, cntL, Rotation, Translation, dy=0):
	z = Translation[2]
	r = radius

	beg = cntP

	# ~~~~~~~~ MIDDLE ~~~~~~~~
	points=[]
	line=[]
	for dt in range(param.dc+1):
		point = curve_point(param, dt, radius, Rotation)
		points.append(point + Translation)
		line.append(cntP)
		cntP+=1
	geo.All_points.append(points)
	cntL = add_c_line(geo, param, line, cntL, dy)

	last = cntP

	ref_point = geo.All_points[-1][0].copy()

	# ~~~~~~~~ LEFT ~~~~~~~~
	points=[]
	line=[]

	line.append(last-1)
	for dl in range(1,param.dl+1):
		point = limb_point_end(param, dl, radius, Rotation)
		points.append(point + Translation)
		line.append(cntP)
		cntP+=1

	geo.All_points.append(points)
	cntL = add_x_line(geo, param, line, cntL, dy, 1, curve_limb=1)

	# ~~~~~~~~ RIGHT ~~~~~~~~
	points=[]
	line=[]

	for dl in range(param.dl):
		point = limb_point_beg(param, dl, radius, Rotation)
		if param.Lshape==0:
			points.append(point + Translation)
		else:
			points.append(point + ref_point)
		line.append(cntP)
		cntP+=1

	line.append(beg)
	geo.All_points.append(points)
	cntL = add_x_line(geo, param, line, cntL, dy, 1, curve_limb=1)

	return cntP, cntL
