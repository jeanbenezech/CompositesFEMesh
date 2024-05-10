import numpy as np
import glob

from utils.parameters import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modified to be called as a library, with optional inputs
def write_inp(E11 = 115.6, E22 = 9.24, nu12 = 0.335, nu23 = 0.487, G12 = 4.826, t_ply = 0.196, K = 1.0e10, K_rig = 1e6, x_spring = 27.5, x_spring_fix = [], x_spring_load = [], K_ground = [], load = -10.0313, displacement = -4.0, rotation_offset = 0.0, StackSeq = [], init_inc = 1.0, min_inc = 1.0e-5, max_inc = 1.0, apply_load = True, fix_to_ground = True):

	param = parameters()
	param.init('parameters')
	nb_ply=param.nbp
	name=param.AbaqusName
	
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
	if not fix_to_ground:
		f.write('100000, {}, 75.0, {}\n'.format(x_spring_fix, -rotation_offset)) # Ground node (need for reaction force)
	f.write('100001, {}, 75.0, 0.0\n'.format(x_spring_fix))
	f.write('100002, {}, 75.0, 420.0\n'.format(x_spring_load))
	# Reference points at support
	f.write('100003, {}, 75.0, {}\n'.format(x_spring_fix, -rotation_offset))
	f.write('100004, {}, 75.0, {}\n'.format(x_spring_load, 420.0+rotation_offset))
	if not fix_to_ground:
		f.write('*NSET, nset=Ground\n')
	f.write('100000\n')
	f.write('*NSET, nset=Fixed_RP\n')
	f.write('100001,\n')
	f.write('*NSET, nset=Load_RP\n')
	f.write('100002,\n')
	f.write('*NSET, nset=Fixed_Sup_RP\n')
	f.write('100003,\n')
	f.write('*NSET, nset=Load_Sup_RP\n')
	f.write('100004,\n')
	# Also need to define spring elements and properties
	f.write('*ELEMENT, type=Spring1, ELSET=Fixed_Spring-spring\n')
	f.write('300001, 100003\n')
	f.write('*ELEMENT, type=Spring1, ELSET=Load_Spring-spring\n')
	f.write('300002, 100004\n')
	# Might need to move to same section as material assignment?
	f.write('*SPRING, elset=Fixed_Spring-spring\n')
	f.write('5\n')
	f.write('{}\n'.format(K))
	f.write('*SPRING, elset=Load_Spring-spring\n')
	f.write('5\n')
	f.write('{}\n'.format(K))
	# Spring for tying base to ground
	if not fix_to_ground:
		f.write('*ELEMENT, type=Spring2, ELSET=Ground_spring\n')
		f.write('300003, 100000, 100003\n')
		f.write('*SPRING, elset=Ground_spring\n')
		f.write('3, 3\n')
		f.write('{}\n'.format(K_ground))

	# Create connector elements
	f.write('*ELEMENT, type=CONN3D2, ELSET=fixed_conn\n')
	f.write('100001, 100001, 100003\n')
	f.write('*ELEMENT, type=CONN3D2, ELSET=load_conn\n')
	f.write('100002, 100002, 100004\n')
	f.write('**\n')

	f.write('*ORIENTATION, name=fixed_end_coord\n')
	f.write('0., 0., -1., 0., 1., 0.\n')
	f.write('1, 0.\n')
	f.write('*ORIENTATION, name=load_end_coord\n')
	f.write('0., 0., 1., 0., 1., 0.\n')
	f.write('1, 0.\n')
	f.write('**\n')

	# Define local  coordinate system
		# ~~~~~~ MATERIAL ASSIGNEMENT ~~~~~~
	# Csection
	# OFFSET = 0 places the reference surface at the centre of the shell,
	# set to -0.5 for the bottom, 0.5 for the top. 
	
	# Note that I have swapped the 45 and -45 degree plies around, as Jean's coordinate 
	# system is orthogonal than that used to lay up the spar
	# Write shell section parameter

	if len(StackSeq) < 1:
		StackSeq = [45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0]

	f.write('*SHELL SECTION, ELSET=Hexahedron, COMPOSITE, ORIENTATION=ori_loc, OFFSET=0, LAYUP="Inner Skin", SYMMETRIC\n')
	for i, ply in enumerate(StackSeq):
		# Subtract 90 for consistency with non-shell version
		f.write('{}, 3, IM7-8552, {}, Ply{}\n'.format(t_ply, ply+90.0, i))

	# Create connector sections
	f.write('**\n')
	f.write('*CONNECTOR SECTION, elset=fixed_conn, behavior=Supp_Trans\n')
	f.write('TRANSLATOR,\n')
	f.write("fixed_end_coord\n",)
	f.write('*CONNECTOR SECTION, elset=load_conn, behavior=Supp_Trans\n')
	f.write('TRANSLATOR,\n')
	f.write("load_end_coord\n")
	f.write('**\n')
	f.write('*end part\n')

	# ~~~~~~~~~~~~~ MATERIAL ~~~~~~~~~~~~~
	# Csection
	f.write('**\n')
	f.write('**---------- MATERIALS ---------- \n')
	f.write('*MATERIAL, NAME=IM7-8552 \n')
	f.write('*ELASTIC, TYPE=ENGINEERING CONSTANTS \n')
	# Give units in GPa for equivalence with kN and mm (for N and mm would be in MPa)
	E33 = E22
	nu13 = nu12
	G13 = G12
	G23 = E22/(2*(1+nu23))
	f.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(E11, E22, E33, nu12, nu13, nu23, G12, G13))
	f.write('{}\n'.format(G23))

	# Define connector behaviour
	f.write('*CONNECTOR BEHAVIOR, name=Supp_Trans\n')
	f.write('*Connector Constitutive Reference\n')
	f.write('{}, , , , , \n'.format(rotation_offset))
	f.write('*Connector Elasticity, component=1\n')
	f.write('{},\n'.format(K_rig))

	# ~~~~~~~~~~~~~ ASSEMBLY ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- ASSEMBLY ---------- \n')
	f.write('*assembly, name=assembly\n')
	f.write('*instance, name=m, part=main\n')
	f.write('*end instance\n')
	f.write('**\n')

	# ~~~~~~~~~ MPCS and Springs ~~~~~~~~~
	# Csection
	# Create New nodes at which above which rotational springs are to be defined
	# Define MPCs
	f.write('*MPC\n')
	f.write('BEAM, m.NSET4, m.Fixed_RP\n')
	f.write('*MPC\n')
	f.write('BEAM, m.NSET5, m.Load_RP\n')	
	f.write('*end assembly\n')

	# ~~~~~~~~~~~~~ INITIAL BOUNDARIES ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- INIT BOUND ---------- \n')
	f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	f.write('m.Fixed_Sup_RP, 1, 1\n')
	f.write('m.Fixed_Sup_RP, 2, 2\n')
	if fix_to_ground:
		f.write('m.Fixed_Sup_RP, 3, 3\n')
	else:
		f.write('m.Ground, 3, 3\n')
	f.write('m.Fixed_Sup_RP, 4, 4\n')
	f.write('m.Fixed_Sup_RP, 6, 6\n')
	f.write('m.Load_Sup_RP, 1, 1\n')
	f.write('m.Load_Sup_RP, 2, 2\n')
	f.write('m.Load_Sup_RP, 4, 4\n')
	f.write('m.Load_Sup_RP, 6, 6\n')
	# Extra ones to turn on/off to restrict twist at the ends
	# f.write('m.Fixed_RP, 4, 4\n') # We might want to allow rotation about the x axis to account for misalignment in the y direction though...
	# f.write('m.Fixed_RP, 6, 6\n')
	# f.write('m.load_RP, 4, 4\n') # We might want to allow rotation about the x axis to account for misalignment in the y direction though...
	#f.write('m.load_RP, 6, 6\n')

	# ~~~~~~~~~~~~~ STEP ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- STEP ---------- \n')
	f.write('*STEP, NAME=Step-1, NLGEOM=YES, SOLVER=ITERATIVE, inc=100\n')
	f.write('*STATIC\n')
	# f.write('0.1, 1.0, 1e-03, 1.0\n')
	f.write('{}, 1.0, {}, {}\n'.format(init_inc, min_inc, max_inc))
	#f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	# Csection
	if apply_load:
		# Apply as point load?
		f.write('*CLOAD \n')
		# Units in kN consistent with above modulus definintion in kN/mm^2 (GPa)
		f.write('m.Load_Sup_RP, 3, {}\n'.format(load))
	else:
		# Otherwise apply as displacement
		f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
		f.write('m.Load_Sup_RP, 3, 3, {}\n'.format(displacement))

	# ~~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- OUTPUT ---------- \n')
	f.write('*OUTPUT, FIELD\n')
	f.write('*NODE OUTPUT\n')
	f.write('U, RF\n')
	f.write('*ELEMENT OUTPUT, elset=m.All_elements\n')
	f.write('E, S\n')

	f.write('*OUTPUT, HISTORY\n')
	f.write('*End Step\n')

if __name__ == '__main__':
	write_inp()