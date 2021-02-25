import numpy as np
import copy

from utils.parameters import *
from utils.geometry import *
from utils.coords import *
from utils.gmsh_hierarchy import *

# ~~~~~~~~~~~~~~ ADD POINTS ~~~~~~~~~~~~~~

def straight_points(geo, param, cntP, cntL, rep, Ylines, stype, dy=0):
	points=[]
	line=[]

	if stype==0:
		point = copy.deepcopy(rep.T)
		points.append(point)
		line.append(cntP)
		cntP+=1
	else:
		line.append(cntP-1)

	# rep.T=points[-1]
	param.substr[stype].Ori = rep.R

	point = translation_along_line(param.substr[stype])
	points.append(point + rep.T)
	line.append(cntP)
	Ylines[1][stype+1] = cntP
	cntP+=1

	rep.T = copy.deepcopy(points[-1])

	geo.All_points.append(points)
	cntL = add_x_line(geo, param, line, cntL, dy)

	return cntP, cntL


def curved_points(geo, param, cntP, cntL, rep, Ylines, stype, savecenter=0):
	points=[]
	line=[]

	if stype==0:
		point = copy.deepcopy(rep.T)
		points.append(point)
		line.append(cntP)
		cntP+=1
	else:
		line.append(cntP-1)

	for dt in range(1, param.dc+1):
		point, center = curve_point(param, dt, rep, stype)
		points.append(point + rep.T)
		line.append(cntP)
		Ylines[1][stype+1] = cntP
		cntP+=1

	if savecenter:
		param.refOriPoint.append(center+rep.T)
	
	rep.T = copy.deepcopy(points[-1])
	rep.R += param.substr[stype].Ori

	geo.All_points.append(points)
	cntL = add_c_line(geo, param, line, cntL, stype)

	return cntP, cntL
