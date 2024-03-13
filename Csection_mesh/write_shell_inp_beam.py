import numpy as np
import glob

from utils.parameters import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modified to be called as a library, with optional inputs
def write_inp(E11 = 115.6, E22 = 9.24, nu12 = 0.335, nu23 = 0.487, G12 = 4.826, t_ply = 0.196, K = 1.0e10, E_steel = 210.0, x_spring = 27.5, x_spring_fix = [], x_spring_load = [], load = -10.0313, displacement = -4.0, rotation_offset = 0.0, StackSeq = [], init_inc = 1.0, min_inc = 1.0e-5, max_inc = 1.0, apply_load = True):

	param = parameters()
	param.init('parameters')
	nb_ply=param.nbp
	name=param.AbaqusName

	# Nodes on fixed end of the spar (maybe pass as input?)
	fixed_end_nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 25, 26, 27, 28, 29, 30, 31, 32,
					33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52,
					53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
					73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
					93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112,
					113, 114, 115, 116, 117, 118, 119, 120, 121, 122] 
	
	# Nodes on loaded end of the spar
	load_end_nodes = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 123, 124, 125, 126, 127, 128, 129, 130, 
				   131, 132, 133, 134, 135, 136, 137, 138, 471, 472, 473, 474, 475, 476, 643, 644, 645, 646, 647, 648,
				   649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668,
				   669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688,
				   689, 690, 691, 692, 693, 694, 695, 696, 863, 864, 865, 866, 867, 868, 1035, 1036, 1037, 1038, 1039, 1040,
				   1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050]

	f=open('Abaqus/'+name+'.inp','w')

	# ~~~~~~~~~~~~~ HEADER ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- INIT ---------- \n')
	f.write('*HEADING\n')
	f.write('Abaqus input file\n')
	f.write('*DISTRIBUTION TABLE, NAME=ori_tab\n')
	f.write('coord3d, coord3d\n')

	# ~~~~~~~~~~~~~ PART ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- PARTS ---------- \n')
	f.write('*part, name=main\n')
	f.write('*INCLUDE, INPUT=Abaqus\\'+param.name+'_mesh.inp\n')
	f.write('*INCLUDE, INPUT=Abaqus\\'+param.name+'_ori.inp\n')
	# Now define reference points and beam elements
	# Define Reference points 
	f.write('*NODE\n')
	if (not x_spring_fix) and (not x_spring_load):
		print("using symmetric offest")
		x_spring_fix = x_spring
		x_spring_load = x_spring
	# Reference points at spar ends
	f.write('100001, {}, 75.0, 0.0\n'.format(x_spring_fix))
	f.write('100002, {}, 75.0, 420.0\n'.format(x_spring_load))
	# Reference points at support
	f.write('100003, {}, 75.0, {}\n'.format(x_spring_fix, -rotation_offset))
	f.write('100004, {}, 75.0, {}\n'.format(x_spring_load, 420.0+rotation_offset))
	f.write('*NSET, nset=Fixed_RP\n')
	f.write('100001,\n')
	f.write('*NSET, nset=Load_RP\n')
	f.write('100002,\n')
	f.write('*NSET, nset=Fixed_Sup_RP\n')
	f.write('100003,\n')
	f.write('*NSET, nset=Load_Sup_RP\n')
	f.write('100004,\n')
	f.write('*ELEMENT, type=T3D2, ELSET=Bars\n')
	for i, node in enumerate(fixed_end_nodes):
		f.write('{}, 100003, {}\n'.format(i+100001, node))
	for i, node in enumerate(load_end_nodes):
		f.write('{}, 100004, {}\n'.format(i+200001, node))
	# Also need to define spring elements and properties
	f.write('*ELEMENT, type=Spring1, ELSET=Fixed_Spring-spring\n')
	f.write('300001, 100003\n')
	f.write('*ELEMENT, type=Spring1, elset=Load_Spring-spring\n')
	f.write('300002, 100004\n')
	# Might need to move to same section as material assignment?
	f.write('*SPRING, elset=Fixed_Spring-spring\n')
	f.write('5\n')
	f.write('{}\n'.format(K))
	f.write('*SPRING, elset=Load_Spring-spring\n')
	f.write('5\n')
	f.write('{}\n'.format(K))
	f.write('**\n')

		# ~~~~~~ MATERIAL ASSIGNEMENT ~~~~~~
	# Csection
	# OFFSET = 0 places the reference surface at the centre of the shell,
	# set to -0.5 for the bottom, 0.5 for the top. 
	# Apparently the specified offset is ignored when using continuum shell eleements
	# For continuum shell elements the ply thickness is only used to assess the fraction
	# of the overall thickness corresponding to each ply. The thickness is taken from the
	# nodal coordinates
	# Continuum shell elements are specified in the mesh definition, which is in a separate
	# .inp file
	
	# Note that I have swapped the 45 and -45 degree plies around, as Jean's coordinate 
	# system is orthogonal than that used to lay up the spar
	if len(StackSeq) < 1:
		StackSeq = [45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0]

	f.write('*SHELL SECTION, ELSET=All_elements, COMPOSITE, ORIENTATION=ori_loc , OFFSET=0, LAYUP="Inner Skin", SYMMETRIC\n')
	for i, ply in enumerate(StackSeq):
		# Subtract 90 for consistency with non-shell version
		f.write('{}, 3, AS4-8552, {}, Ply{}\n'.format(t_ply, ply+90.0, i))

	# Write shell section parameter
	f.write('*SOLID SECTION, ELSET=Bars, MATERIAL=Steel\n')
	f.write('{}\n'.format(29400.0/len(fixed_end_nodes))) # Cross-sectional area (use CSA of end blocks in mm^2 divided by number of bars)
	f.write('**\n')
	f.write('*end part\n')

	# ~~~~~~~~~~~~~ MATERIAL ~~~~~~~~~~~~~
	# Csection
	f.write('**\n')
	f.write('**---------- MATERIALS ---------- \n')
	f.write('*MATERIAL, NAME=AS4-8552 \n')
	f.write('*ELASTIC, TYPE=ENGINEERING CONSTANTS \n')
	# Give units in GPa for equivalence with kN and mm (for N and mm would be in MPa)
	E33 = E22
	nu13 = nu12
	G13 = G12
	G23 = E22/(2*(1+nu23))
	f.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(E11, E22, E33, nu12, nu13, nu23, G12, G13))
	f.write('{}\n'.format(G23))
	f.write('*MATERIAL, NAME=Steel \n')
	f.write('*ELASTIC, TYPE=ISOTROPIC \n')
	f.write('{}, 0.3\n'.format(E_steel)) # Modulus, Poisson's ratio

	# ~~~~~~~~~~~~~ ASSEMBLY ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- ASSEMBLY ---------- \n')
	f.write('*assembly, name=assembly\n')
	f.write('*instance, name=m, part=main\n')
	f.write('*end instance\n')
	f.write('**\n')

	# ~~~~~~~~~ MPCS and Springs ~~~~~~~~~
	# Csection
	
	# Used if MPC is applied as Kinematic Coupling constraint
	# f.write('*SURFACE, TYPE=NODE, NAME=jean_nset4, internal\n')
	# f.write('m.nset4, 1.\n')
	# f.write('*COUPLING, CONSTRAINT NAME=BotSurfKC, REF NODE=m.MasterNode4, SURFACE=jean_nset4\n')
	# f.write('*KINEMATIC\n')
	# f.write('*SURFACE, TYPE=NODE, NAME=jean_nset5, internal\n')
	# f.write('m.nset5, 1.\n')
	# f.write('*COUPLING, CONSTRAINT NAME=TopSurfKC, REF NODE=m.MasterNode5, SURFACE=jean_nset5\n')
	#f.write('*KINEMATIC\n')
	
	# Create New nodes at which above which rotational springs are to be defined
	# Define MPCs
	f.write('*MPC\n')
	f.write('BEAM, m.NSET4, m.Fixed_RP\n')
	f.write('*MPC\n')
	f.write('BEAM, m.NSET5, m.Load_RP\n')

	# Write MPCs at lead applicaton point uing an equation constraint
	#f.write('*EQUATION\n')
	#f.write('2\n')
	#f.write('m.nset5, 3, 1.0, m.MasterNode5, 3, -1.0\n')
	
	f.write('*end assembly\n')

	# ~~~~~~~~~~~~~ INITIAL BOUNDARIES ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- INIT BOUND ---------- \n')
	f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	f.write('Fixed_RP, 1, 1\n')
	f.write('Fixed_RP, 2, 2\n')
	f.write('Fixed_RP, 3, 3\n')
	f.write('Fixed_RP, 4, 4\n')
	f.write('Fixed_RP, 6, 6\n')
	f.write('Load_RP, 1, 1\n')
	f.write('Load_RP, 2, 2\n')
	f.write('Load_RP, 4, 4\n')
	f.write('Load_RP, 6, 6\n')
	
	# Csection
	# f.write('m.MasterNode4, 1, 1, 0.0\n')
	# f.write('m.MasterNode4, 2, 2, 0.0\n')
	# f.write('m.MasterNode4, 3, 3, 0.0\n')
	# f.write('m.MasterNode5, 1, 1, 0.0\n')
	# f.write('m.MasterNode5, 2, 2, 0.0\n')
	# Laminate

	# ~~~~~~~~~~~~~ STEP ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- STEP ---------- \n')
	f.write('*STEP, NAME=Step-1, NLGEOM=YES, SOLVER=ITERATIVE, inc=100\n')
	f.write('*STATIC\n')
	# f.write('0.1, 1.0, 1e-03, 1.0\n')
	f.write('{}, 1.0, {}, {}\n'.format(init_inc, min_inc, max_inc))
	#f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	# Csection
	#f.write('m.MasterNode5, 3, 3, -1\n')

	f.write('*CLOAD \n')
	# Units in kN consistent with above modulus definintion in kN/mm^2 (GPa)
	f.write('Load_RP, 3, {}\n'.format(load))

	# ~~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- OUTPUT ---------- \n')
	f.write('*OUTPUT, FIELD\n')
	f.write('*NODE OUTPUT\n')
	f.write('U, RF\n')
	f.write('*ELEMENT OUTPUT, elset=m.All_elements\n')
	f.write('E, S\n')

	f.write('*OUTPUT, HISTORY\n')
	# The below lines change the default parameters associated with incremenation control for assessing convergence.
	# These settings will typically allow more increments than the default and so I think are there to force 
	# f.write('*CONTROLS, PARAMETERS=FIELD\n')
#	f.write('5e-2,,,,,,,\n')
	# f.write('5,5,,,,,,\n')
	# f.write('*CONTROLS, PARAMETERS=TIME INCREMENTATION\n')
#	f.write('1000,1000,1000,1000,1000, , , , ,10,1000\n')
	# f.write('100,100,100,100,100, , , , ,10,100\n')
	# f.write('*CONTROLS, PARAMETERS=LINE SEARCH\n')
	# f.write('10,,,,\n')

	f.write('*End Step\n')

if __name__ == '__main__':
	write_inp()