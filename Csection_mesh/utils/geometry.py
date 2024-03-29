import numpy as np
import copy

from utils.parameters import *

class reference():
	def __init__(self,
				 O=None, # origin
				 T=None, # translation
				 R=None): # rotation
		self.O = np.asarray([0.0, 0.0, 0.0])
		self.T = np.asarray([0.0, 0.0, 0.0])
		self.R = 0 * np.pi/180.
	
	def update_line_ref(self, r=0):
		# self.O = np.asarray([r*np.sin(self.R), -r*np.cos(self.R), 0.0])
		self.O = np.asarray([0, -r, 0.0])
		self.T = copy.deepcopy(self.O)
		self.R = 0 * np.pi/180.


class geometry():
	def __init__(self,
				 All_points=None,
				 Lines=None,
				 Transfinite_curves=None,
				 Surfaces=None,
				 Physical_surfaces=None,
				 Volumes=None,
				 Physical_volumes=None,
				 delta_p=None,
				 delta_l=None,
				 delta_s=None,
				 delta_v=None,
				 go_to_side=None):
		self.All_points=All_points
		self.Lines=Lines
		self.Surfaces=Surfaces

	def init(self, param):
		self.delta_p=0
		self.delta_l=0
		self.delta_s=0
		self.delta_v=0

		self.go_to_side=(param.dy-1)

		self.All_points=[]
		self.Lines=[]
		self.Surfaces=[]
		self.Transfinite_curves=[]

		for i in range(2): # 1 vector for the number of the line and one for the discretization of the line
			self.Transfinite_curves.append([])

		# Surfaces groups for boundary conditions
		nb_physical_surf=6 # right, left, front, back, bottom and top for the stiffener

		self.Physical_surfaces=[]
		for i in range(nb_physical_surf):
			self.Physical_surfaces.append([])

		# Volumes groups for the stacking sequence and the matrix domain
		self.Physical_volumes=[]

		len_PG=param.nbp

		param.resin_from_StackSeq = 0
		for v in param.StackSeq:
			if v[2] > 0:
				param.resin_from_StackSeq = max(param.resin_from_StackSeq,int(v[2]))
		len_PG += param.resin_from_StackSeq  # +1 for each the resin layer group

		if param.resin==1:
			len_PG+=1 # +1 for the resin layer group
		if param.CZ==1:
			len_PG+=1 # +1 for the cohesize elements

		for i in range(len_PG):
			self.Physical_volumes.append([])
		# print(len(self.Physical_volumes))

	# def ajust(self, param): # TODO fill this function if needed
		# for i in range(len(self.All_points)):
		# 	for j in range(len(self.All_points[i])):
		# 		if param.type_==2:
		# 			self.All_points[i][j][1]+=param.X
		# 			self.All_points[i][j][0]+=param.limb
		# 		else:
		# 			self.All_points[i][j][1]+=param.Y
		# 			if param.type_==1:
		# 				self.All_points[i][j][0]+=param.limb
		# 				if param.flat_limb:
		# 					self.All_points[i][j][0]+=param.side

	def dump(self, param):
		f=open(param.name+'.geo', 'w')
		f.write('lc = '+str(param.lc)+';\n')

		cnt=0
		# Points
		# print(self.All_points)
		for i in range(len(self.All_points)):
			for j in range(len(self.All_points[i])):
				f.write('Point('+str(cnt)+') = {'+str(self.All_points[i][j][0])+', ' \
				+str(self.All_points[i][j][1])+', '+str(self.All_points[i][j][2])+', lc};\n') # write each points
				cnt+=1

		# self.Lines
		for i in range(len(self.Lines)):
			if len(self.Lines[i])>1:
				#~ f.write('Line('+str(i+1)+') = {'+str(self.Lines[i][0]))
				f.write('Spline('+str(i+1)+') = {'+str(self.Lines[i][0]))
				for j in range(1,len(self.Lines[i])):
					f.write(', '+str(self.Lines[i][j]))
				f.write('};\n')

		# Transfinite curves
		for i, l in enumerate(self.Transfinite_curves[0]):
			f.write('Transfinite Curve {'+str(l+1)+'} = '+str(self.Transfinite_curves[1][i])+';\n')


		cnt_empty_surf = 0
		# Surfaces
		for i in range(len(self.Surfaces)):
			if len(self.Surfaces[i]) > 1:
				f.write('Curve Loop('+str(i+1)+') = {'+str(self.Surfaces[i][0]+1)+', '+ \
				str(self.Surfaces[i][1]+1)+', -'+str(self.Surfaces[i][3]+1)+', -'+str(self.Surfaces[i][2]+1)+'};\n')
				f.write('Surface('+str(i+1)+') = {'+str(i+1)+'};\n')
			else:
				cnt_empty_surf+=1

		# Transfinite Surfaces
		f.write('Transfinite Surface {'+str(1))
		for i in range(1, len(self.Surfaces)):
			if len(self.Surfaces[i]) > 1:
				f.write(', '+str(i+1))
		f.write('};\n')

		if param.recombine==1:
			# Recombine Surfaces --> To mesh quads
			f.write('Recombine Surface {'+str(1))
			for i in range(1, len(self.Surfaces)):
				if len(self.Surfaces[i]) > 1:
					f.write(', '+str(i+1))
			f.write('};\n')

		# ~~~~~~~~~ EXTRUSION ~~~~~~~~~
		f.write('Extrude {0, 0, '+str(param.Z)+'} {Surface{'+str(1)+'}')
		for i in range(1, len(self.Surfaces)):
			if len(self.Surfaces[i]) > 1:
				f.write('; Surface{'+str(i+1)+'}')
		if len(param.dz)==1:
			f.write('; Layers{'+str(param.dz[0])+'}; Recombine; }\n')
		elif param.ramp==0:
			f.write('; Layers{ {'+str(param.dz[0]))
			for dz in param.dz[1:]:
				f.write(', '+str(dz))
			f.write('}, {'+str(1/len(param.dz)))
			i=1
			for dz in param.dz[1:-1]:
				f.write(', '+str((1+i)/len(param.dz)))
				i+=1
			f.write(', 1} }; Recombine; }\n')
		else: # Ramp + len(dz)>1
			assert len(param.dz) == len(param.dz_intervals)
			f.write('; Layers{ {'+str(param.dz[0]))
			for dz in param.dz[1:]:
				f.write(', '+str(dz))
			f.write('}, {'+str(param.dz_intervals[0]))
			for interval in param.dz_intervals[1:]:
				f.write(', '+str(interval))
			f.write('} }; Recombine; }\n')

		# # ~~~~~~~~~ PHYSICAL GROUPS ~~~~~~~~~

		# delta_pg = 0

		# Physical volumes
		# print (self.Physical_volumes)
		for i, j in enumerate(self.Physical_volumes):
			if len(j)>0:
				f.write('Physical Volume('+str(i+1)+') = {'+str(j[0]+1))
				for k in j[1:]:
					f.write(', '+str(k+1))
				f.write('};\n')


		# # print (self.Physical_surfaces)
		# # Physical surfaces
		# for i, j in enumerate(self.Physical_surfaces):
		# 	if len(j)>0:
		# 		f.write('Physical Surface('+str(i+1+delta_pg)+') = {'+str(j[0]+1))
		# 		for k in j[1:]:
		# 			f.write(', '+str(k+1))
		# 		f.write('};\n')
