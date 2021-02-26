import numpy as np

from utils.parameter import *
from utils.geometry import *
from utils.points import *
from utils.gmsh_hierarchy import *

def flat_line_beg(param, geo, Translation):

	cntP=geo.delta_p
	cntL=geo.delta_l
	cntS=geo.delta_s

	y = -param.Y					# y-coord to begin with
	delta_Y = param.Y/(param.nbp) # ply thickness

	ptype=0

	last_comp_layer = param.dy-1
	if param.CZ==1:
		last_comp_layer = param.dy-2

	# First line
	cntP, cntL = Flat_points(geo, param, cntP, cntL, Translation, y)
	if(param.flat_limb):
		cntP, cntL = Flat_sides(geo, param, cntP, cntL, Translation, y)
		# Zlines which not exist but to keep Lines iterator
		for z_incr in range(6):
			geo.Lines.append([0])
			cntL+=1
	else:
		for z_incr in range(4):
			geo.Lines.append([0])
			cntL+=1

	for dy in range(1, param.dy):
		if param.resin==0:
			y += delta_Y # comp layer thickness
		else:
			if param.CZ==0:
				if dy%2==1: # comp layer thickness
					y += delta_Y-param.e
					ptype=0
				else: # resin thickness
					y += param.e
					ptype=1
			else:
				rT = param.e
				cT = param.e / 100.
				pT = delta_Y - 2*rT - cT

				if dy%4==1: # comp layer thickness
					y += pT
					ptype=0
				elif dy%2==0: # resin thickness
					y += rT
					ptype=1
				elif dy%4==3: # cohezive thickness
					y += cT
					ptype=2

		if(param.flat_limb==0):
			cntL = add_y_lines(geo, param, cntL, dy)
		else:

			if dy<last_comp_layer or param.resin==0:
				cntL = add_y_lines_with_side(geo, param, cntL, dy)
			else:
				cntL = add_y_lines_with_side(geo, param, cntL, dy, last=1)

		cntP, cntL = Flat_points(geo, param, cntP, cntL, Translation, y, dy)
		if(param.flat_limb):
			if param.resin==1:
				if dy<param.dy-1: # For CZ, the last segment (param.dy-2) is add just for the construction but no surfaces are link to it
					cntP, cntL = Flat_sides(geo, param, cntP, cntL, Translation, y, dy)
				else:
					for z_incr in range(2):
						geo.Lines.append([0])
						cntL+=1
			else:
				cntP, cntL = Flat_sides(geo, param, cntP, cntL, Translation, y, dy)
			# Zlines which not exist but to keep Lines iterator
			for z_incr in range(6):
				geo.Lines.append([0])
				cntL+=1
		else:
			for z_incr in range(4):
				geo.Lines.append([0])
				cntL+=1

		if(param.flat_limb==0):
			cntS = front_surf(param, geo, cntL, cntS, dy, layer_type=ptype)
		else:
			if dy<last_comp_layer or param.resin==0:
				cntS = front_surf_with_side(param, geo, cntL, cntS, dy, layer_type=ptype)
			else:
				cntS = front_surf_with_side(param, geo, cntL, cntS, dy, layer_type=ptype, last=1)

	geo.delta_p=cntP
	geo.delta_l=cntL
	geo.delta_s=cntS
	geo.delta_v+=(param.nbp)

	if(param.flat_limb==0):
		return geo.Lines[-7][0],geo.Lines[-7][-1], geo.Lines[-6][-1], geo.Lines[-5][0], geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v
	else:
		return geo.Lines[-11][0],geo.Lines[-11][-1], geo.Lines[-10][-1], geo.Lines[-9][0], geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v

def curved_line_beg(param, geo, Translation, Rotation, record_volume=1):
	cntP=geo.delta_p
	cntL=geo.delta_l
	cntS=geo.delta_s

	max_radius = param.X
	min_radius = param.X-param.R

	ptype=0

	radius = min_radius
	delta_R = param.R/(param.nbp_c) # Ply thickness

	# First line
	radius = min_radius
	cntP, cntL = Lshape_points(geo, param, radius, cntP, cntL, Rotation, Translation, 0)

	for z_incr in range(4):
		geo.Lines.append([0])
		cntL+=1

	for dy in range(1, param.dy_c):

		if param.resin==0:
			radius += delta_R # comp layer thickness
		else:
			if param.CZ==0:
				if dy%2==1:
					radius += delta_R-param.e # comp layer thickness
					ptype=0
				else:
					radius += param.e  # resin thickness
					ptype=1
			else:
				rT = param.e
				cT = param.e / 100.
				pT = delta_R - 2*rT - cT

				if dy%4==1: # comp layer thickness
					radius += pT
					ptype=0
				elif dy%2==0: # resin thickness
					radius += rT
					ptype=1
				elif dy%4==3: # cohezive thickness
					radius += cT
					ptype=2


		cntL = add_y_lines(geo, param, cntL, dy, bool_curve=1)

		cntP, cntL = Lshape_points(geo, param, radius, cntP, cntL, Rotation, Translation, dy)

		for z_incr in range(4):
			geo.Lines.append([0])
			cntL+=1

		if record_volume==1:
			cntS = front_surf(param, geo, cntL, cntS, dy, bool_curve=1, layer_type=ptype)
		# else:
		# 	geo.Surfaces.append([])
		# 	geo.Surfaces.append([])
		# 	geo.Surfaces.append([])
		# 	cntS+=3

	geo.delta_p=cntP
	geo.delta_l=cntL
	geo.delta_s=cntS
	geo.delta_v+=(param.nbp_c)

	return geo.Lines[-7][0],geo.Lines[-7][-1], geo.Lines[-6][-1], geo.Lines[-5][0], geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v
