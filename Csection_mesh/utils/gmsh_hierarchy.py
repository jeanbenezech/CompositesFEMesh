import numpy as np

from utils.parameters import *
from utils.geometry import *

# ~~~~~~~~~~~~~~ ADD LINES ~~~~~~~~~~~~~~

def add_x_line(geo, param, line, cntL, dy=0):
	geo.Lines.append(line) #0 -8 ;;; #5 -3
	geo.Transfinite_curves[0].append(cntL)
	geo.Transfinite_curves[1].append(param.dx) # TODO discretization function of the substr

	cntL+=1
	return cntL

def add_c_line(geo, param, line, cntL, stype):
	geo.Lines.append(line) #0 -8 ;;; #5 -3
	geo.Transfinite_curves[0].append(cntL)
	val = int(param.dc)#* (1.-(90.-param.substr[stype].Ori)/90.))
	geo.Transfinite_curves[1].append(val)

	cntL+=1
	return cntL

def add_y_lines(geo, param, cntL, Ylines, layer_type=0):

	if layer_type==0: # ply layer
		val=param.ddy
	elif layer_type==1: # resin layer
		val=param.ddy
	elif layer_type==2: # Cohezive layer
		val=1

	for i in range(len(param.substr)+1):
		geo.Lines.append([Ylines[0][i], Ylines[1][i]])
		geo.Transfinite_curves[0].append(cntL)
		geo.Transfinite_curves[1].append(val)
		cntL+=1

	return cntL

# ~~~~~~~~~~~~~~ ADD SURFACES ~~~~~~~~~~~~~~

def front_surf(param, geo, rep, cntL, cntS, dy, layer_type=0):

	incr_vol = int((dy-1)/2)
	if param.CZ and param.resin:
		incr_vol=int((dy-1)/4)
	# print(geo.delta_v)
	if param.resin==0 and param.CZ==0: # Case where there is no resin layer between plies
		incr_vol = (dy-1)

	start = cntL - (4*len(param.substr) + 2)

	for k in range(0,len(param.substr)):
		yline_start = 3* len(param.substr)
		newXline_start = 2* len(param.substr)
		geo.Surfaces.append([start + k, start + k + yline_start + 2, start + k + yline_start + 1, start + k + newXline_start + 1])

		if layer_type==1: # resin layer
			if param.CZ:
				geo.Physical_volumes[-2].append(cntS)
			else:
				geo.Physical_volumes[-1].append(cntS)
		elif layer_type==2: # Cohezive layer
			geo.Physical_volumes[-1].append(cntS)
		else:
			geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS)

		cntS+=1

	return cntS