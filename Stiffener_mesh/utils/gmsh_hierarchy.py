import numpy as np

from utils.parameter import *
from utils.geometry import *

# ~~~~~~~~~~~~~~ ADD LINES ~~~~~~~~~~~~~~

def add_x_line(geo, param, line, cntL, dy=0, add_in_filler=0, position=1, curve_limb=0, bool_side=0):
	geo.Lines.append(line) #0 -8 ;;; #5 -3
	geo.Transfinite_curves[0].append(cntL)
	if position==0:
		geo.Transfinite_curves[1].append(param.dx)
	else:
		if bool_side==0:
			geo.Transfinite_curves[1].append(param.dl)
		else:
			geo.Transfinite_curves[1].append(param.dside)


	if dy == param.dy-1 and add_in_filler and curve_limb==0:
		geo.Filler_edge[0].append(cntL)

	if curve_limb==1 and add_in_filler and dy == param.dy_c-1:
		geo.Filler_edge[0].append(cntL)

	cntL+=1
	return cntL

def add_c_line(geo, param, line, cntL, dy=0):
	geo.Lines.append(line) #0 -8 ;;; #5 -3
	geo.Transfinite_curves[0].append(cntL)
	val = int(param.dc * (1.-(90.-param.theta)/90.))
	geo.Transfinite_curves[1].append(val)

	if dy == param.dy_c-1:
		geo.Filler_edge[0].append(cntL)

	cntL+=1
	return cntL

def add_y_lines(geo, param, cntL, dy, bool_curve=0):

	if param.resin==0:
		val=param.ddy
	else:
		if param.CZ==1:
			if dy%4==1: # comp layer discr
				val=param.ddy
			elif dy%2==0:  # resin layer discr
				val=param.ddy_r
			else:
				val=1
		else:
			if dy%2==1:
				val=param.ddy
			else:
				val=param.ddy_r

	Npbl  = param.dx+1
	Npb_limb  = param.dl
	Npb3l = param.dx+1 + param.dl*2

	if bool_curve:
		Npbl  = param.dc+1
		Npb_limb  = param.dl
		Npb3l = param.dc+1 + param.dl*2


	# ~~~~~ Ligne 7 ~~~~~
	geo.Lines.append([(dy-1)*Npb3l + Npbl-1 + geo.delta_p, dy*Npb3l + Npbl-1 + geo.delta_p]) #7
	geo.Transfinite_curves[0].append(cntL)
	geo.Transfinite_curves[1].append(val)

	# ~~~~~ Ligne 8 ~~~~~
	geo.Lines.append([(dy-1)*Npb3l + geo.delta_p, dy*Npb3l + geo.delta_p]) #8
	geo.Transfinite_curves[0].append(cntL+1)
	geo.Transfinite_curves[1].append(val)

	# if bool_curve==0:
	# ~~~~~ Ligne 9 ~~~~~
	geo.Lines.append([(dy-1)*Npb3l + Npb_limb+Npbl-1 + geo.delta_p, dy*Npb3l + Npb_limb+Npbl-1 + geo.delta_p]) #9
	geo.Transfinite_curves[0].append(cntL+2)
	geo.Transfinite_curves[1].append(val)

	# ~~~~~ Ligne 10 ~~~~~
	geo.Lines.append([(dy-1)*Npb3l + Npb_limb+Npbl-1 + geo.delta_p +1, dy*Npb3l + Npb_limb+Npbl-1 + geo.delta_p +1]) #10
	geo.Transfinite_curves[0].append(cntL+3)
	geo.Transfinite_curves[1].append(val)
	# else:
	# 	# ~~~~~ Ligne 9 ~~~~~
	# 	geo.Lines.append([(dy-1)*Npb3l + 2*Npb_limb+Npbl-1 + geo.delta_p, dy*Npb3l + 2*Npb_limb+Npbl-1 + geo.delta_p]) #9
	# 	geo.Transfinite_curves[0].append(cntL+2)
	# 	geo.Transfinite_curves[1].append(val)
	#
	# 	# ~~~~~ Ligne 10 ~~~~~
	# 	geo.Lines.append([(dy-1)*Npb3l + Npbl-1 + geo.delta_p +1, dy*Npb3l + Npbl-1 + geo.delta_p +1]) #10
	# 	geo.Transfinite_curves[0].append(cntL+3)
	# 	geo.Transfinite_curves[1].append(val)

	cntL+=4
	return cntL

def add_y_lines_with_side(geo, param, cntL, dy, last=0):

	Npbl  = param.dx+1
	Npb_limb  = param.dl
	Npb_side  = param.dl
	Npb5l = param.dx+1 + param.dl*2 + param.dl*2

	if param.resin==0:
		val=param.ddy
	else:
		if param.CZ==1:
			if dy%4==1: # comp layer discr
				val=param.ddy
			elif dy%2==0:  # resin layer discr
				val=param.ddy_r
			else:
				val=1
		else:
			if dy%2==1:
				val=param.ddy
			else:
				val=param.ddy_r

	# ~~~~~ Ligne 11 ~~~~~
	geo.Lines.append([(dy-1)*Npb5l + Npbl-1 + geo.delta_p, dy*Npb5l + Npbl-1 + geo.delta_p]) #7
	geo.Transfinite_curves[0].append(cntL)
	geo.Transfinite_curves[1].append(val)

	# ~~~~~ Ligne 12 ~~~~~
	geo.Lines.append([(dy-1)*Npb5l + geo.delta_p, dy*Npb5l + geo.delta_p]) #8
	geo.Transfinite_curves[0].append(cntL+1)
	geo.Transfinite_curves[1].append(val)

	# ~~~~~ Ligne 13 ~~~~~
	geo.Lines.append([(dy-1)*Npb5l + (Npb_side+Npbl-1) + geo.delta_p, dy*Npb5l + Npb_side+Npbl-1 + geo.delta_p]) #9
	geo.Transfinite_curves[0].append(cntL+2)
	geo.Transfinite_curves[1].append(val)

	# ~~~~~ Ligne 14 ~~~~~
	geo.Lines.append([(dy-1)*Npb5l + Npb_side+Npbl-1 + geo.delta_p +1, dy*Npb5l + Npb_side+Npbl-1 + geo.delta_p +1]) #10
	geo.Transfinite_curves[0].append(cntL+3)
	geo.Transfinite_curves[1].append(val)

	if last==0:
		# ~~~~~ Ligne 15 ~~~~~
		geo.Lines.append([(dy-1)*Npb5l + (2*Npb_side+Npb_limb+Npbl-1) + geo.delta_p, dy*Npb5l + 2*Npb_side+Npb_limb+Npbl-1 + geo.delta_p]) #9
		geo.Transfinite_curves[0].append(cntL+4)
		geo.Transfinite_curves[1].append(val)

		# ~~~~~ Ligne 16 ~~~~~
		geo.Lines.append([(dy-1)*Npb5l + 2*Npb_side+Npb_limb+Npbl-1 + geo.delta_p +1, dy*Npb5l + 2*Npb_side+Npb_limb+Npbl-1 + geo.delta_p +1]) #10
		geo.Transfinite_curves[0].append(cntL+5)
		geo.Transfinite_curves[1].append(val)
	else:
		geo.Lines.append([0])
		geo.Lines.append([0])

	cntL+=6

	return cntL

# ~~~~~~~~~~~~~~ ADD SURFACES ~~~~~~~~~~~~~~

def front_surf(param, geo, cntL, cntS, dy, bool_curve=0, layer_type=0):

	incr_vol = int((dy-1)/2)
	if param.CZ:
		incr_vol=int((dy-1)/4)
	# print(geo.delta_v)
	if param.resin==0: # Case where there is no resin layer between plies
		incr_vol = (dy-1)

	cmpt0 = cntL-18 # 18 = nb line by surface

	geo.Surfaces.append([cmpt0, cmpt0+7, cmpt0+8, cmpt0+11])
	geo.Surfaces.append([cmpt0+1, cmpt0+9, cmpt0+7, cmpt0+12])
	geo.Surfaces.append([cmpt0+2, cmpt0+8, cmpt0+10, cmpt0+13])

	if layer_type==1: # resin layer
		geo.Physical_volumes[-2].append(cntS)
		geo.Physical_volumes[-2].append(cntS+1)
		geo.Physical_volumes[-2].append(cntS+2)
	elif layer_type==2: # Cohezive layer
		geo.Physical_volumes[-3].append(cntS)
		geo.Physical_volumes[-3].append(cntS+1)
		geo.Physical_volumes[-3].append(cntS+2)
	else:
		geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS)
		geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS+1)
		geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS+2)

	cntS+=3

	return cntS

def front_surf_with_side(param, geo, cntL, cntS, dy, layer_type=0, last=0):

	incr_vol = int((dy-1)/2)
	if param.CZ:
		incr_vol=int((dy-1)/4)
	if param.resin==0:
		incr_vol = (dy-1)

	cmpt0 = cntL-28 # 28 = nb line by surface

	geo.Surfaces.append([cmpt0, cmpt0+11, cmpt0+12, cmpt0+17])
	geo.Surfaces.append([cmpt0+1, cmpt0+13, cmpt0+11, cmpt0+18])
	geo.Surfaces.append([cmpt0+2, cmpt0+12, cmpt0+14, cmpt0+19])
	if last==0:
		geo.Surfaces.append([cmpt0+3, cmpt0+15, cmpt0+13, cmpt0+20])
		geo.Surfaces.append([cmpt0+4, cmpt0+14, cmpt0+16, cmpt0+21])

	if layer_type==1: # resin layer
		geo.Physical_volumes[-2].append(cntS)
		geo.Physical_volumes[-2].append(cntS+1)
		geo.Physical_volumes[-2].append(cntS+2)
		if last==0:
			geo.Physical_volumes[-2].append(cntS+3)
			geo.Physical_volumes[-2].append(cntS+4)
	elif layer_type==2: # Cohezive layer
		geo.Physical_volumes[-3].append(cntS)
		geo.Physical_volumes[-3].append(cntS+1)
		geo.Physical_volumes[-3].append(cntS+2)
		if last==0:
			geo.Physical_volumes[-3].append(cntS+3)
			geo.Physical_volumes[-3].append(cntS+4)
	else:
		geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS)
		geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS+1)
		geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS+2)
		geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS+3)
		geo.Physical_volumes[incr_vol+geo.delta_v].append(cntS+4)

	cntS+=5
	if last==1:
		cntS-=2
	return cntS

def front_inter_surf(geo, param, cntL, cntS, dy):

	incr_vol = dy-1

	move_to_next_surf = 3 * (move_to_next_surf)

	geo.Surfaces.append([cntL-8, cntL-5, cntL-4, cntL-3]) # (5)  [delta+0, delta+3, delta+4, delta+5]
	geo.Volumes[incr_vol+geo.delta_v].append(cntS)

	geo.Volumes[incr_vol+geo.delta_v+move_to_next_surf].append(cntS)

	cntS+=1

	return cntS
