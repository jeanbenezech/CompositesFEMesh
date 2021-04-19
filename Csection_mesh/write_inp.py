import numpy as np
import glob

from utils.parameters import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

	param = parameters()
	param.init()

	fix=[0]
	pull=1
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
	for i in range(1,nb_ply+1):
		f.write('*SOLID SECTION, ELSET=elset'+str(i)+', ORIENTATION=ori_loc, MATERIAL=D-8552\n')

	f.write('*COHESIVE SECTION, ELSET=elset'+str(nb_ply+1)+', MATERIAL=CZ, RESPONSE=TRACTION SEPARATION, THICKNESS=GEOMETRY \n')
	# ~~~~~~ KINEMATIC COUPLING ~~~~~~

	for fi in fix:
		f.write('*KINEMATIC COUPLING, REF NODE=MasterNode'+str(fi)+'\n')
		f.write('nset'+str(fi)+', 1, 6\n')
	f.write('*KINEMATIC COUPLING, REF NODE=MasterNode'+str(pull)+'\n')
	f.write('nset'+str(pull)+', 1, 6\n')
	f.write('*end part\n')


	# ~~~~~~~~~~~~~ MATERIAL ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- MATERIALS ---------- \n')
	f.write('*MATERIAL, NAME=D-8552 \n')
	f.write('*ELASTIC, TYPE=ENGINEERING CONSTANTS \n')
	f.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(137300., 8800., 8800., 0.314, 0.314, 0.487, 4900., 4900.))
	f.write('{}\n'.format(2960.))
	#~ f.write('*Expansion, type=ORTHO \n')
	#~ f.write('-2.1e-07, 3.3e-05, 3.3e-05 \n')
	# f.write('*MATERIAL, name=FILLER \n')
	# f.write('*ELASTIC, type=ENGINEERING CONSTANTS \n')
	# f.write('8800.,  8800., 137300.,  0.487,   0.02,   0.02,  2960.,  4900. \n')
	# f.write('4900., \n')
	#~ f.write('*Expansion, type=ORTHO \n')
	#~ f.write('3.3e-05,  3.3e-05, -2.1e-07 \n')
	# f.write('*MATERIAL, name=PEEK \n')
	# f.write('*ELASTIC \n')
	# f.write('2500., 0.4 \n')
	#~ f.write('*Expansion \n')
	#~ f.write('4.7e-05, \n')

	f.write('*MATERIAL, name=CZ \n')
	f.write('*ELASTIC, TYPE=TRACTION \n')
	f.write(' 7612.,1370.,1370.\n')
	f.write('*DAMAGE INITIATION, CRITERION=QUADS \n')
	f.write(' 74.2, 110.4, 110.4 \n')
	f.write('*DAMAGE EVOLUTION, TYPE=ENERGY, MIXED MODE BEHAVIOR=BK, POWER=1.45 \n')
	f.write(' 0.3, 0.87, 0.87 \n')
	# f.write('*DAMAGE STABILIZATION \n')
	# f.write('0.00001 \n')

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
	for fi in fix:
		f.write('m.MasterNode'+str(fi)+', 1, 1, 0.0\n')
		f.write('m.MasterNode'+str(fi)+', 2, 2, 0.0\n')
		f.write('m.MasterNode'+str(fi)+', 3, 3, 0.0\n')


	# ~~~~~~~~~~~~~ STEP ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- STEP ---------- \n')
	f.write('*STEP, NAME=Step-1, SOLVER=ITERATIVE, inc=5000\n')
	f.write('*STATIC, stabilize, factor=1e-2\n')
	# f.write('0.1, 1.0, 1e-03, 1.0\n')
	f.write('1e-02, 1.0, 1e-06, 5e-01\n')
	f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	f.write('m.MasterNode'+str(pull)+', 2, 2, -2\n')
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
