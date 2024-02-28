import numpy as np
import glob

from utils.parameters import *

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

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
	f.write('**\n')

	# ~~~~~~~~~~~~~ PART ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- PARTS ---------- \n')
	f.write('*part, name=main\n')
	f.write('*INCLUDE, INPUT='+param.name+'_mesh.inp\n')
	f.write('*INCLUDE, INPUT='+param.name+'_ori.inp\n')
	# ~~~~~~ MATERIAL ASSIGNEMENT ~~~~~~
	# Laminate
	f.write('*SOLID SECTION, ELSET=elset1, ORIENTATION=ori_loc, MATERIAL=UD\n')
	f.write('*SOLID SECTION, ELSET=elset3, ORIENTATION=ori_loc, MATERIAL=UD\n')
	f.write('*SOLID SECTION, ELSET=elset4, ORIENTATION=ori_loc, MATERIAL=UD\n')
	f.write('*SOLID SECTION, ELSET=elset5, ORIENTATION=ori_loc, MATERIAL=RESIN\n')
	f.write('*SOLID SECTION, ELSET=elset6, ORIENTATION=ori_loc, MATERIAL=DELAM\n')

	# ~~~~~~ KINEMATIC COUPLING ~~~~~~
	# Laminate
	f.write('*KINEMATIC COUPLING, REF NODE=MasterNode4\n') # ZMIN clamped
	f.write('nset4, 1, 6\n')
	f.write('*KINEMATIC COUPLING, REF NODE=MasterNode2\n')  # YMIN dirichlet
	f.write('nset2, 1, 6\n')
	f.write('*end part\n')


	# ~~~~~~~~~~~~~ MATERIAL ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- MATERIALS ---------- \n')
	f.write('*MATERIAL, NAME=UD \n')
	f.write('*ELASTIC, TYPE=ENGINEERING CONSTANTS \n')
	f.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(
		137100., 8800., 8800., 0.314, 0.314, 0.487, 4900., 4900.))
	f.write('{}\n'.format(2959.))
	f.write('*MATERIAL, NAME=RESIN \n')
	f.write('*ELASTIC, TYPE=ISOTROPIC \n')
	f.write('{}, {}\n'.format(8800., 0.314))
	f.write('*MATERIAL, NAME=DELAM \n')
	f.write('*ELASTIC, TYPE=ISOTROPIC \n')
	f.write('{}, {}\n'.format(0.01, 0.314))


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
	f.write('m.MasterNode2, 1, 1, 0.0\n')
	f.write('m.MasterNode2, 2, 2, 0.0\n')
	# Laminate


	# ~~~~~~~~~~~~~ STEP ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- STEP ---------- \n')
	f.write('*STEP, NAME=Step-1, SOLVER=ITERATIVE, inc=5000\n')
	f.write('*STATIC, stabilize, factor=1e-2\n')
	f.write('1e-02, 1.0, 1e-06, 5e-01\n')
	f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	f.write('m.MasterNode2, 3, 3, 10\n')
	# ~~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- OUTPUT ---------- \n')
	f.write('*OUTPUT, FIELD\n')
	f.write('*NODE OUTPUT\n')
	f.write('U, RF\n')
	f.write('*ELEMENT OUTPUT, elset=m.All_elements\n')
	f.write('EVOL, E, S\n')

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
