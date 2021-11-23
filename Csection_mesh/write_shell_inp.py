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
	# Csection
	# MUST BE CHANGED
	f.write('*SHELL SECTION, ELSET=All_elements, COMPOSITE, ORIENTATION=ori_loc , OFFSET=0, LAYUP="Inner Skin", SYMMETRIC\n')
	f.write('0.20833333333333334, 3, AS4-8552, 45., Ply1\n')
	f.write('0.20833333333333334, 3, AS4-8552, -45., Ply2\n')
	f.write('0.20833333333333334, 3, AS4-8552, 45., Ply3\n')
	f.write('0.20833333333333334, 3, AS4-8552, -45., Ply4\n')
	f.write('0.20833333333333334, 3, AS4-8552, 45., Ply5\n')
	f.write('0.20833333333333334, 3, AS4-8552, -45., Ply6\n')
	f.write('0.20833333333333334, 3, AS4-8552, 0., Ply7\n')
	f.write('0.20833333333333334, 3, AS4-8552, 90., Ply8\n')
	f.write('0.20833333333333334, 3, AS4-8552, 0., Ply9\n')
	f.write('0.20833333333333334, 3, AS4-8552, 90., Ply10\n')
	f.write('0.20833333333333334, 3, AS4-8552, 0., Ply11\n')
	f.write('0.20833333333333334, 3, AS4-8552, 90., Ply12\n')
	f.write('*end part\n')

	# ~~~~~~~~~~~~~ MATERIAL ~~~~~~~~~~~~~
	# Csection
	f.write('**\n')
	f.write('**---------- MATERIALS ---------- \n')
	f.write('*MATERIAL, NAME=AS4-8552 \n')
	f.write('*ELASTIC, TYPE=ENGINEERING CONSTANTS \n')
	f.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(137300., 8800., 8800., 0.314, 0.314, 0.487, 4900., 4900.))
	f.write('{}\n'.format(2960.))


	# ~~~~~~~~~~~~~ ASSEMBLY ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- ASSEMBLY ---------- \n')
	f.write('*assembly, name=assembly\n')
	f.write('*instance, name=m, part=main\n')
	f.write('*end instance\n')
	# ~~~~~~ KINEMATIC COUPLING ~~~~~~
	# Csection
	f.write('*SURFACE, TYPE=NODE, NAME=jean_nset4, internal\n')
	f.write('m.nset4, 1.\n')
	f.write('*COUPLING, CONSTRAINT NAME=BotSurfKC, REF NODE=m.MasterNode4, SURFACE=jean_nset4\n')
	f.write('*KINEMATIC\n')
	f.write('*SURFACE, TYPE=NODE, NAME=jean_nset5, internal\n')
	f.write('m.nset5, 1.\n')
	f.write('*COUPLING, CONSTRAINT NAME=TopSurfKC, REF NODE=m.MasterNode5, SURFACE=jean_nset5\n')
	f.write('*KINEMATIC\n')
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
	f.write('*STEP, NAME=Step-1, NLGEOM=YES, SOLVER=ITERATIVE, inc=5000\n')
	f.write('*STATIC\n')
	# f.write('0.1, 1.0, 1e-03, 1.0\n')
	f.write('1e-01, 1.0, 1e-06, 5e-01\n')
	f.write('*BOUNDARY, TYPE=DISPLACEMENT\n')
	# Csection
	f.write('m.MasterNode5, 3, 3, -1\n')
	# Laminate
	# ~~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~~
	f.write('**\n')
	f.write('**---------- OUTPUT ---------- \n')
	f.write('*OUTPUT, FIELD\n')
	f.write('*NODE OUTPUT\n')
	f.write('U, RF\n')
	f.write('*ELEMENT OUTPUT, elset=m.All_elements\n')
	f.write('E, S\n')

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
