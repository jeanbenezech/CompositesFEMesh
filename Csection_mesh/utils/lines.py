import numpy as np
import copy

from utils.parameters import *
from utils.geometry import *
from utils.points import *
from utils.gmsh_hierarchy import *

def Clines(param, geo, rep):
	cntP=geo.delta_p
	cntL=geo.delta_l
	cntS=geo.delta_s

	ptype=0

	# delta_R = param.R/(param.nbp) # Ply thickness


	Ylines=[np.zeros(len(param.substr)+1, dtype=int),np.zeros(len(param.substr)+1, dtype=int)]

	rep.update_line_ref(0)
	Ylines[1][0] = cntP
	for stype in range(0,len(param.substr)):
		if param.substr[stype].shape=='straight':
			cntP, cntL = straight_points(geo, param, cntP, cntL, rep, Ylines, stype)
		elif param.substr[stype].shape=='curved':
			cntP, cntL = curved_points(geo, param, cntP, cntL, rep, Ylines, stype, savecenter=1)


	for incr in range(len(param.substr)+1): # TODO ensure the lines compt for surface creation : the first line have no yline linking it to the previous line, as no previous exist
		geo.Lines.append([0])
		cntL+=1
	Ylines[0] = copy.deepcopy(Ylines[1])


	radius = 0
	for dy in range(1, param.dy):

		if param.resin==1 or param.CZ ==1:
			delta_R = param.StackSeq[int((dy-1)/2)][1]


		if param.resin==0 and param.CZ == 0:
			radius += param.StackSeq[(dy-1)][1] # comp layer thickness
			ptype = param.StackSeq[(dy-1)][2]

		elif param.resin==1 and param.CZ == 0:
			if dy%2==1:
				radius += delta_R-param.e # comp layer thickness
				ptype=0
			else:
				radius += param.e  # resin thickness
				ptype=1

		elif param.resin==0 and param.CZ ==1:
			cT = param.e / 10.
			if dy%2==1:
				radius += delta_R-cT # comp layer thickness
				ptype=0
			else:
				radius += cT  # cohezive thickness
				ptype=2

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

		rep.update_line_ref(radius)
		Ylines[1][0] = cntP

		for stype in range(0,len(param.substr)):
			param.substr[stype].thick = radius
			if param.substr[stype].shape=='straight':
				cntP, cntL = straight_points(geo, param, cntP, cntL, rep, Ylines, stype)
			elif param.substr[stype].shape=='curved':
				cntP, cntL = curved_points(geo, param, cntP, cntL, rep, Ylines, stype)
		
		cntL = add_y_lines(geo, param, cntL, Ylines, layer_type=ptype)
		Ylines[0] = copy.deepcopy(Ylines[1])

		cntS = front_surf(param, geo, rep, cntL, cntS, dy, layer_type=ptype)

	geo.delta_p=cntP
	geo.delta_l=cntL
	geo.delta_s=cntS
	geo.delta_v+=(param.nbp)

	return geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v
