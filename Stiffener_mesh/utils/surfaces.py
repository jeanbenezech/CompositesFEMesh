import numpy as np

from utils.geometry import *
from utils.parameter import *
from utils.lines import *

def front_surface(param, geo):

	N=np.zeros(12)


	# Base
	Translation = np.asarray([0.0, 0.0, 0.0])
	N[0], N[1], N[6], N[7], geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v = flat_line_beg(param, geo, Translation)

	# Curve left
	# if param.resin:
		# Translation = np.asarray([0.0, param.X-param.e, 0.0])
	Translation = np.asarray([0.0, param.X+param.fl, 0.0])

	# else:
	# 	Translation = np.asarray([0.0, param.X+param.fl, 0.0])

	if param.Lshape==1:
		limit_angle=90

	Rotation = 0
	N[2], N[3], N[8], N[9], geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v = curved_line_beg(param, geo, Translation, Rotation)

	if param.Lshape==0:
		record_volume=1
		# Curve right
		# if param.resin:
			# Translation = np.asarray([2.0*(param.X), param.X-param.e, 0.0])
		Translation = np.asarray([2.0*(param.X)+param.fl, param.X+param.fl, 0.0])
		# else:
		# 	Translation = np.asarray([2.0*(param.X+param.e), param.X+param.fl, 0.0])
		# Translation = np.asarray([2.0*param.X, param.X-param.e, 0.0])
		Rotation = np.pi/2
		N[4], N[5], N[10], N[11], geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v = curved_line_beg(param, geo, Translation, Rotation, record_volume)

	return geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v, N

def filler(param, geo, N):
	# Filler boundaries FRONT
	geo.Lines.append([N[0], N[3]])
	geo.Lines.append([N[2], N[5]])
	geo.Lines.append([N[4], N[1]])

	# Filler boundaries FRONT EXT
	geo.Lines.append([N[6], N[11]])
	geo.Lines.append([N[8], N[7]])
	geo.Lines.append([N[10], N[9]])

	# Faire en sorte qu'il n'y ait pas plusieurs éléments dans l'epaisseur des zones de resine
	# entre la flat plate et les L-shapes
	for i in range(6):
		geo.Transfinite_curves[0].append(geo.delta_l+i)
		geo.Transfinite_curves[1].append(param.ddy_r)

	geo.Filler_surf=[]
	# les +1 et -1 sont la pour la numerotation finale des noeuds : besoin d'avoir les noeuds entre 1 et n (et pas 0 n-1 comme le vecteur Lines)

	# Filler
	geo.Filler_surf.append([ geo.Filler_edge[0][0]+1, -(len(geo.Lines)-6)-1,  geo.Filler_edge[0][3]+1, -(len(geo.Lines)-5)-1,  geo.Filler_edge[0][6]+1, -(len(geo.Lines)-4)-1])

	# print(geo.Filler_surf)
	# print(geo.Filler_edge)

	s=[]
	for i in range(12):
		s.append([])

	# print(len(geo.Lines)) # =299
	# s3 = {295, 288, 299, 216};

	# s6 = {296, 138, 297, 289};

	# s9 = {294, 215, 298, 139};


	s[3].append((len(geo.Lines)-5)+1)  # 295

	s[6].append((len(geo.Lines)-4)+1)  # 296

	s[9].append((len(geo.Lines)-6)+1)  # 294

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	s[3].append( geo.Filler_edge[0][1+3*2]+1) # 288

	s[6].append( geo.Filler_edge[0][1+3*0]+1) # 138

	s[9].append( geo.Filler_edge[0][1+3*1]+1) # 215

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	s[3].append((len(geo.Lines)-1)+1)  # 299

	s[6].append((len(geo.Lines)-3)+1)  # 297

	s[9].append((len(geo.Lines)-2)+1)  # 298

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	s[3].append( geo.Filler_edge[0][2+3*1]+1) # 216

	s[6].append( geo.Filler_edge[0][2+3*2]+1) # 289

	s[9].append( geo.Filler_edge[0][2+3*0]+1) # 139

	# for i in range(12):
	# print(s[3])
	geo.Filler_surf.append(s[3])
	geo.Filler_surf.append(s[6])
	geo.Filler_surf.append(s[9])

def filler_L_shape(param, geo, N):

	cntL = geo.delta_l
	cntP = geo.delta_p

	# pointN5 = [geo.,y,0.0]
	flat_All_points=[]

	for i in range(len(geo.All_points)):
		for j in range(len(geo.All_points[i])):
			flat_All_points.append(geo.All_points[i][j].copy())

	N5_point = flat_All_points[int(N[2])]
	N5_point[0] += param.fl

	N4_point = flat_All_points[int(N[1])]
	N4_point[1] += param.fl

	center = []
	center.append(N5_point[0])
	center.append(N4_point[1])
	points=[]
	line=[]
	x1 = N4_point[0] - center[0]
	y1 = N5_point[1] - center[1]
	for dt in range(param.dc+1):
		lame = 1./2.
		x = np.power(dt/param.dc, 1./lame) * x1
		y = y1/x1 * np.power(np.power(x1, lame) - np.power(x, lame), 1./lame)
		point = [center[0] + x, center[1] + y, 0.0]
		points.append(point)
		line.append(cntP)
		cntP+=1
	geo.All_points.append(points)

	geo.Lines.append(line)
	geo.Transfinite_curves[0].append(cntL)
	geo.Transfinite_curves[1].append(param.dc)
	cntL+=1


	# Filler boundaries FRONT
	geo.Lines.append([N[0], N[3]])
	val = geo.Lines[-1-1][0]
	geo.Lines.append([N[2], val])
	val = geo.Lines[-1-2][-1]
	geo.Lines.append([val, N[1]])

	# Filler boundaries FRONT EXT
	# geo.Lines.append([N[6], N[11]])
	geo.Lines.append([N[8], N[7]])
	# geo.Lines.append([N[10], N[9]])

	# Faire en sorte qu'il n'y ait pas plusieurs éléments dans l'epaisseur des zones de resine
	# entre la flat plate et les L-shapes
	for i in range(4):
		geo.Transfinite_curves[0].append(cntL+i)
		geo.Transfinite_curves[1].append(param.ddy_r)

	geo.Filler_surf=[]
	# les +1 et -1 sont la pour la numerotation finale des noeuds : besoin d'avoir les noeuds entre 1 et n (et pas 0 n-1 comme le vecteur Lines)

	# Filler
	geo.Filler_surf.append([ geo.Filler_edge[0][0]+1, -(len(geo.Lines)-4)-1,  geo.Filler_edge[0][3]+1, -(len(geo.Lines)-3)-1,  -(len(geo.Lines)-5)-1, -(len(geo.Lines)-2)-1])

	# print(geo.Filler_surf)
	# print(geo.Filler_edge)

	s=[]

	# print(len(geo.Lines)) # =299
	# s3 = {295, 288, 299, 216};

	# s6 = {296, 138, 297, 289};

	# s9 = {294, 215, 298, 139};

	s.append((len(geo.Lines)-4)+1)  # 294

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	s.append( geo.Filler_edge[0][1+3*1]+1) # 215

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	s.append((len(geo.Lines)-1)+1)  # 298

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	s.append( geo.Filler_edge[0][2+3*0]+1) # 139


	geo.Filler_surf.append(s)
