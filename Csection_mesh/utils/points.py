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

	newRep = rep.T
	# rep.T=points[-1]
	param.substr[stype].Ori = rep.R

	if stype==0:
		### Flat rotated flange
		if(param.substr[stype].type=='flangeRotRi' or param.substr[stype].type=='flangeRotLe'):
			Ix = param.substr[stype].local_Ix0(param, -rep.O[1])
			param.substr[stype].L = Ix
			### modify rep to have a 0 flat rotated surf
			diff = param.substr[stype].L
			newRep[0] += diff

		point = copy.deepcopy(newRep)
		points.append(point)
		line.append(cntP)
		cntP+=1

		### Flat rotated flange
		if(param.substr[stype].type=='flangeRotRi' or param.substr[stype].type=='flangeRotLe'):
			Ix = param.substr[stype].local_Ix0(param, -rep.O[1])
			param.substr[stype].L = Ix
			### modify rep to have a 0 flat rotated surf
			diff = param.substr[stype].L
			newRep[0] -= diff
	else:
		line.append(cntP-1)



	if(param.substr[stype].type=='flangeRotRi' or param.substr[stype].type=='flangeRotLe'):
		if(param.substr[stype].type=='flangeRotLe'):
			Ix = param.substr[stype].local_Ix(param, -rep.O[1])
			# if stype>0:
			param.substr[stype].L = Ix - param.substr[stype-1].A
			# else:
			# 	param.substr[stype].L = Ix
			# 	### modify rep to have a 0 flat rotated surf
			# 	diff = param.substr[stype].L - param.substr[stype].Linit
			# 	newRep[0] += diff
		if(param.substr[stype].type=='flangeRotRi'):
			Ix = param.substr[stype].local_Ix(param, -rep.O[1])
			if stype==len(param.dflange_intervals)+4: ## +4 for the substruc of the top of the spar
				param.substr[stype].L = param.Height - Ix ## first substurcture going down
			else:
				param.substr[stype].L = param.substr[stype-1].A - Ix

	point = translation_along_line(param.substr[stype])
	points.append(point + newRep)
	line.append(cntP)
	Ylines[1][stype+1] = cntP
	cntP+=1

	rep.T = copy.deepcopy(points[-1])

	if(param.substr[stype].type=='flangeRotRi' or param.substr[stype].type=='flangeRotLe'):
		# print(stype)
		# print(param.substr[stype].FlangeRotAngle)
		# local_y = param.substr[stype].local_y(param.R-rep.O[1])
		if(param.substr[stype].type=='flangeRotLe'):
			Ix = param.substr[stype].local_Ix(param, -rep.O[1])
			# if stype>0:
			param.substr[stype].L = Ix - param.substr[stype-1].A
			diff = param.substr[stype].L - param.substr[stype].Linit
			rep.T[0] -= diff
			# else:
			# 	param.substr[stype].L = Ix
			# 	diff = param.substr[stype].L - param.substr[stype].Linit
			# 	rep.T[0] -= 2.*diff
		if(param.substr[stype].type=='flangeRotRi'):
			Ix = param.substr[stype].local_Ix(param, -rep.O[1])
			if stype==len(param.dflange_intervals)+4: ## +4 for the substruc of the top of the spar
				param.substr[stype].L = param.Height - Ix ## first substurcture going down
			elif stype>len(param.dflange_intervals)+4:
				param.substr[stype].L = param.substr[stype-1].A - Ix
			diff = param.substr[stype].L - param.substr[stype].Linit
			rep.T[0] += diff

	geo.All_points.append(points)
	cntL = add_x_line(geo, param, line, cntL, stype)

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
