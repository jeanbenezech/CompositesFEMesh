import numpy as np
import glob

from utils.parameters import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

	param = parameters()
	param.init()

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
	f.write('**\n')

	# ~~~~~~~~~~~~~ PART ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- PARTS ---------- \n')
	f.write('*part, name=main\n')
	f.write('*INCLUDE, INPUT='+param.name+'_mesh.inp\n')
	f.write('*INCLUDE, INPUT='+param.name+'_ori.inp\n')
	# ~~~~~~ MATERIAL ASSIGNEMENT ~~~~~~
	# Csection
	f.write('*SOLID SECTION, ELSET=Hexahedra, ORIENTATION=ori_loc, MATERIAL=AS4-8552\n')
	f.write('*COHESIVE SECTION, ELSET=Cohesive, MATERIAL=CZ, RESPONSE=TRACTION SEPARATION, THICKNESS=GEOMETRY \n')
	# Laminate
	# for i in range(1,nb_ply+1):
	# 	f.write('*SOLID SECTION, ELSET=elset'+str(i)+', ORIENTATION=ori_loc, MATERIAL=D-8552\n')
	# f.write('*COHESIVE SECTION, ELSET=elset'+str(nb_ply+1)+', MATERIAL=CZ, RESPONSE=TRACTION SEPARATION, THICKNESS=GEOMETRY \n')

	# ~~~~~~ KINEMATIC COUPLING ~~~~~~
	# Csection
	f.write('*KINEMATIC COUPLING, REF NODE=MasterNode4\n')
	f.write('nset4, 1, 6\n')
	f.write('*KINEMATIC COUPLING, REF NODE=MasterNode5\n')
	f.write('nset5, 1, 6\n')
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
	f.write('**---------- MATERIALS ---------- \n')
	f.write('*MATERIAL, NAME=AS4-8552 \n')
	f.write('*ELASTIC, TYPE=ENGINEERING CONSTANTS \n')
	f.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(137300., 8800., 8800., 0.314, 0.314, 0.487, 4900., 4900.))
	f.write('{}\n'.format(2960.))

	# Csection
	f.write('*MATERIAL, name=CZ \n')
	f.write('*ELASTIC, TYPE=TRACTION \n')
	f.write(' 7612.,1370.,1370.\n')
	f.write('*DAMAGE INITIATION, CRITERION=QUADS \n')
	f.write(' 74.2, 110.4, 110.4 \n')
	f.write('*DAMAGE EVOLUTION, TYPE=ENERGY, MIXED MODE BEHAVIOR=BK, POWER=1.45 \n')
	f.write(' 0.3, 0.87, 0.87 \n')

	# ~~~~~~~~~~~~~ ASSEMBLY ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- ASSEMBLY ---------- \n')
	f.write('*assembly, name=assembly\n')
	f.write('*instance, name=m, part=main\n')
	f.write('*end instance\n')
	f.write('*end assembly\n')

	# ~~~~~~~~~~~~~ INITIAL BOUNDARIES ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- INIT BOUND ---------- \n')
	f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	# Csection
	f.write('m.MasterNode4, 1, 1, 0.0\n')
	f.write('m.MasterNode4, 2, 2, 0.0\n')
	f.write('m.MasterNode4, 3, 3, 0.0\n')
	f.write('m.MasterNode5, 1, 1, 0.0\n')
	f.write('m.MasterNode5, 2, 2, 0.0\n')
	# Laminate


	# ~~~~~~~~~~~~~ STEP ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- STEP ---------- \n')
	f.write('*STEP, NAME=Step-1, SOLVER=ITERATIVE, inc=5000\n')
	f.write('*STATIC, stabilize, factor=1e-2\n')
	# f.write('0.1, 1.0, 1e-03, 1.0\n')
	f.write('1e-02, 1.0, 1e-06, 5e-01\n')
	f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	# Csection
	f.write('m.MasterNode5, 3, 3, -2\n')
	# Laminate
	# ~~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- OUTPUT ---------- \n')
	f.write('*OUTPUT, FIELD\n')
	f.write('*NODE OUTPUT\n')
	f.write('U, RF\n')
	f.write('*ELEMENT OUTPUT, elset=m.All_elements\n')
	f.write('EVOL, E, S\n')
	f.write('*ELEMENT OUTPUT, elset=m.Cohesive\n')
	f.write('SDEG, QUADSCRT\n')

	f.write('*OUTPUT, HISTORY\n')
	f.write('*CONTROLS, PARAMETERS=FIELD\n')
#	f.write('5e-2,,,,,,,\n')
	f.write('5,5,,,,,,\n')
	f.write('*CONTROLS, PARAMETERS=TIME INCREMENTATION\n')
#	f.write('1000,1000,1000,1000,1000, , , , ,10,1000\n')
	f.write('100,100,100,100,100, , , , ,10,100\n')
	f.write('*CONTROLS, PARAMETERS=LINE SEARCH\n')
	f.write('10,,,,\n')

	f.write('*End Step\n')
