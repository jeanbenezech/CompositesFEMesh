import numpy as np
import glob

from utils.parameters import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if __name__ == '__main__':
# Modified to be called as a library, with optional inputs
def write_inp(E11 = 115.6, E22 = 9.24, nu12 = 0.335, nu23 = 0.487, G12 = 4.826, K = 1.0e10, x_spring = 27.5, load = -10.0313, init_inc = 1.0, min_inc = 1.0e-5, max_inc = 1.0):

	param = parameters()
	param.init('parameters')
	nb_ply=param.nbp
	# name=param.AbaqusName From my latest change - keep if it breaks code
	name=param.name

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
	# Assumes Abaqus is being run from the parent directory, and not the Abaqus subfolder
	f.write('*INCLUDE, INPUT=Abaqus\\'+param.name+'_mesh.inp\n')
	f.write('*INCLUDE, INPUT=Abaqus\\'+param.name+'_ori.inp\n')
	# ~~~~~~ MATERIAL ASSIGNEMENT ~~~~~~
	# Layer section
	# f.write('*SOLID SECTION, ELSET=Hexahedra, ORIENTATION=ori_loc, MATERIAL=ALU\n')
	# Cohesive section
	# f.write('*COHESIVE SECTION, ELSET=Cohesive, MATERIAL=CZ, RESPONSE=TRACTION SEPARATION, THICKNESS=GEOMETRY \n')

	# Laminate
	# for i in range(1,nb_ply+1):
	#	f.write('*SOLID SECTION, ELSET=elset'+str(i)+', ORIENTATION=ori_loc, MATERIAL=AS4-8552\n')
	#	f.write('*COHESIVE SECTION, ELSET=elset'+str(nb_ply+1)+', MATERIAL=CZ, RESPONSE=TRACTION SEPARATION, THICKNESS=GEOMETRY \n')
	
	# Alternatively just write single Section applied to all elements
	f.write('*SOLID SECTION, ELSET=Hexahedron, ORIENTATION=ori_loc, MATERIAL=AS4-8552\n')

	# ~~~~~~ KINEMATIC COUPLING ~~~~~~
	# Csection
	f.write('*end part\n')
	# Laminate
	# f.write('*KINEMATIC COUPLING, REF NODE=MasterNode'+str(fi)+'\n')
	# f.write('nset'+str(fi)+', 1, 6\n')
	# f.write('*KINEMATIC COUPLING, REF NODE=MasterNode'+str(pull)+'\n')
	# f.write('nset'+str(pull)+', 1, 6\n')
	# f.write('*end part\n')


	# ~~~~~~~~~~~~~ MATERIAL ~~~~~~~~~~~~~
	# Csection
	f.write('**\n')
	# f.write('**---------- MATERIALS ---------- \n')
	# f.write('*MATERIAL, NAME=ALU \n')
	# f.write('*ELASTIC, TYPE=ISOTROPIC \n')
	# f.write('{}, {}\n'.format(70000., 0.3))

	f.write('*MATERIAL, NAME=AS4-8552 \n')
	f.write('*ELASTIC, TYPE=ENGINEERING CONSTANTS \n')
	# Give units in GPa for equivalence with kN and mm (for N and mm would be in MPa)
	E33 = E22
	nu13 = nu12
	G13 = G12
	G23 = E22/(2*(1+nu23))
	f.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(E11, E22, E33, nu12, nu13, nu23, G12, G13))
	f.write('{}\n'.format(G23))

	# Cohesive behavior
	# f.write('*MATERIAL, name=CZ \n')
	# f.write('*ELASTIC, TYPE=TRACTION \n')
	# # f.write(' 7612.,1370.,1370.\n')
	# f.write(' 2244898.,1250000.,754841.\n')
	# f.write('*DAMAGE INITIATION, CRITERION=QUADS \n')
	# f.write(' 74.2, 110.4, 110.4 \n')
	# f.write('*DAMAGE EVOLUTION, TYPE=ENERGY, MIXED MODE BEHAVIOR=BK, POWER=1.45 \n')
	# f.write(' 0.3, 0.87, 0.87 \n')

	# ~~~~~~~~~~~~~ ASSEMBLY ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- ASSEMBLY ---------- \n')
	f.write('*assembly, name=assembly\n')
	f.write('*instance, name=m, part=main\n')
	f.write('*end instance\n')
	f.write('**\n')

	# ~~~~~~~~~ MPCS and Springs ~~~~~~~~~~
	# Csection

	# Code for writing MPCs as Kinematic coupling constraint
	# f.write('*SURFACE, TYPE=NODE, NAME=surface_nset4, internal\n')
	# f.write('m.nset4, 1.\n')
	# f.write('*SURFACE, TYPE=NODE, NAME=surface_nset5, internal\n')
	# f.write('m.nset5, 1.\n')
	# f.write('*COUPLING, CONSTRAINT NAME=TopSurfKC, REF NODE=m.MasterNode5, SURFACE=surface_nset5\n')
	# f.write('*KINEMATIC\n')
	#f.write('3')

	# Create New nodes at which above which rotational springs are to be defined
	f.write('*Node\n')
	f.write('1, {}, 75.0, 0.0\n'.format(x_spring))
	f.write('2, {}, 75.0, 420.0\n'.format(x_spring))
	f.write('**\n')
	# Define Node sets
	f.write('*Nset, nset=Fixed_RP\n')
	f.write('1,\n')
	f.write('*Nset, nset=Load_RP\n')
	f.write('2,\n')
	# Define MPCs
	f.write('*MPC\n')
	f.write('BEAM, m.NSET4, Fixed_RP\n')
	f.write('*MPC\n')
	f.write('BEAM, m.NSET5, Load_RP\n')
	f.write('*Spring, elset=Fixed_Spring-spring\n')
	f.write('5\n')
	f.write('{}\n'.format(K))
	f.write('*Element, type=Spring1, elset=Fixed_Spring-spring\n')
	f.write('1, 1\n')
	f.write('*Spring, elset=Load_Spring-spring\n')
	f.write('5\n')
	f.write('{}\n'.format(K))
	f.write('*Element, type=Spring1, elset=Load_Spring-spring\n')
	f.write('2, 2\n')

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
	f.write('**\n')
	# Csection
	
	#f.write('m.nset4, 1, 1, 0.0\n')
	#f.write('m.nset4, 2, 2, 0.0\n')
	#f.write('m.nset4, 3, 3, 0.0\n')
	

	#f.write('m.nset5, 1, 1, 0.0\n')
	#f.write('m.nset5, 2, 2, 0.0\n')
	# Laminate


	# ~~~~~~~~~~~~~ STEP ~~~~~~~~~~~~~

	f.write('**---------- STEP ---------- \n')
	f.write('*STEP, NAME=Step-1, NLGEOM=YES, SOLVER=ITERATIVE, inc=100\n')
	f.write('*STATIC\n')
	# f.write('0.1, 1.0, 1e-03, 1.0\n')
	# f.write('1e-02, 1.0, 1e-06, 5e-01\n')
	f.write('{}, 1.0, {}, {}\n'.format(init_inc, min_inc, max_inc))
	# f.write('1.0, 1.0, 1e-05, 1.0\n')
	# Csection
	#f.write('m.nset5, 3, 3, -2\n')

	# f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	f.write('*CLOAD \n')
	# Units in kN consistent with above modulus definintion in kN/mm^2 (GPa)
	f.write('Load_RP, 3, {}\n'.format(load))
	
	# Laminate
	# ~~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- OUTPUT ---------- \n')
	f.write('*OUTPUT, FIELD\n')
	f.write('*NODE OUTPUT\n')
	f.write('U, RF\n')
	f.write('*ELEMENT OUTPUT, elset=m.All_elements\n')
	f.write('EVOL, E, S\n')
	# f.write('*ELEMENT OUTPUT, elset=m.Cohesive\n')
	# f.write('SDEG, QUADSCRT\n')

	f.write('*OUTPUT, HISTORY\n')
	# The below lines change the default parameters associated with incremenation control for assessing convergence.
	# These settings will typically allow more increments than the default and so I think are there to force 
	# convergence. Consider commenting these
	# f.write('*CONTROLS, PARAMETERS=FIELD\n')
	# f.write('5,5,,,,,,\n')
	# f.write('*CONTROLS, PARAMETERS=TIME INCREMENTATION\n')
	# f.write('100,100,100,100,100, , , , ,10,100\n')
	# f.write('*CONTROLS, PARAMETERS=LINE SEARCH\n')
	# f.write('10,,,,\n')

	f.write('*End Step\n')
