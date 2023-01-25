import sys
import random
import numpy as np

# geometry in millimeters

def write_parameters():
	nb_plies = 24
	nb_wrinkles = 0
	Xlenght = 100.0
	Ylenght = 50.0 # thickness
	Zlenght = 100.0
	height = 0.0

	# ply_thickness = Ylenght/(nb_plies+0.0)
	# ply_type = [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0]
	# ply_thickness = [0.15, 0.05, 0.01, 0.05, 0.10, 0.05, 0.01, 0.05, 0.15]
	# ply_thickness = [0.17, 0.03, 0.01, 0.03, 0.14, 0.03, 0.01, 0.03, 0.17]
	# ply_thickness = [0.10, 0.05, 0.025, 0.01, 0.0175, 0.035, 0.07, 0.035, 0.0175, 0.01, 0.025, 0.05, 0.10]
	# tot_ply_thickness = 0.0
	# for t in ply_thickness:
	# 	tot_ply_thickness += t
	# StackSeq = [0.0, 0.0, 0.0, 45.0, 45.0, 45.0, 0.0, -45.0, -45.0]
	# StackSeq = [0.0, 0.0]
	# StackSeq = [0.0, 0.0, 0.0, 0.0, 45.0, 45.0, 45.0, 45.0, 45.0, 0.0, -45.0,-45.0, -45.0]
	StackSeq = [45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 90.0, 0.0, 90.0, 0.0, 90.0, 0.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0]


	tot_ply_thickness = Ylenght
	ply_thickness = np.full_like(StackSeq, Ylenght/(nb_plies+0.0))
	ply_type = np.full_like(StackSeq, 0) # 0 - plies ; 1 - resin // Here: only plies, if resin, it is done via auto resin otherwize need to be specified

	# WRINKLES Parameters
	minWsize = 1.0
	maxWsize = 5.0
	minWposX = 0.30*Xlenght
	maxWposX = 0.70*Xlenght
	minWposY = 0.49* Ylenght
	maxWposY = 0.51* Ylenght
	minWposZ = 0.45*Zlenght
	maxWposZ = 0.55*Zlenght
	minWdampX = 0.01
	maxWdampX = 0.06
	minWdampY = 0.05
	maxWdampY = 0.15
	minWdampZ = 0.05
	maxWdampZ = 0.10
	minWori = -20
	maxWori = 20

	parameters=open('parameters.txt','w')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('name(s)                  : RVE\n') # mesh name
	parameters.write('Shape(i)                 : 1\n')   # 0(default): Csection ; 1: flat laminate
	parameters.write('Auto_Resin_betw_plies(b) : 0\n')   # 1: yes ; 0: no
	parameters.write('cohezive_elements(b)     : 0\n')   # 1: yes ; 0: no
	parameters.write('recombine(b)             : 1\n')   # 1: recombine mesh: hex ;  0: no: only prisms
	parameters.write('Shell(b)                 : 0\n')   # 1: recombine mesh: hex ;  0: no: only prisms
	parameters.write('GaussianThickness(b)     : 0\n')   # 1: gridtrans ;  0: flat
	parameters.write('CornerThickness(b)       : 0\n')   # 1: variation ;  0: flat
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('np(i)         : '+str(nb_plies)+'\n')  #
	parameters.write('X(f)          : '+str(Xlenght)+'\n')   #
	parameters.write('Y(f)          : '+str(tot_ply_thickness)+'\n')  #
	parameters.write('R(f)          : 10\n')    #
	parameters.write('Height(f)     : '+str(height)+'\n')    #
	parameters.write('Z(f)          : '+str(Zlenght)+'\n')   #
	parameters.write('e(f)          : 0.01\n')  # Resin layer thickness
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~STACKSEQ~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	for i in range(nb_plies):
		if i<10:
			parameters.write('p'+str(i)+'(f,f,b)   : '+str(StackSeq[i])+','+str(ply_thickness[i])+','+str(ply_type[i])+'\n')
		else:
			parameters.write('p'+str(i)+'(f,f,b)  : '+str(StackSeq[i])+','+str(ply_thickness[i])+','+str(ply_type[i])+'\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~MESH~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('lc(f)    : 1\n')        # mesh caracteristic size
	parameters.write('dx(i)    : 20\n')       #
	parameters.write('ddy(i)   : 3\n')        # discretisation of a ply through thickness :: 2 = 1 elem/ply 
	parameters.write('dz(i)    : 20\n')       #
	parameters.write('dc(i)    : 0\n')        # corner :: CSPAR
	parameters.write('dflange(i) : 0\n')      # discretization of the flange :: CSPAR
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~TRANSFORMATION~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('wrinkle(i)             : '+str(nb_wrinkles)+'\n') # Do we need to add wrinkle in the gridMod (c++) code, 2+ for multiple wrinkles
	parameters.write('Ramp(b)                : 0\n') # CSPAR
	parameters.write('RotateRVE(b)           : 0\n') #
	parameters.write('AngleRotateRVE(f)      : 45.\n')
	parameters.write('InteriorRadiusRVE(f)   : 100.\n')
	for i in range(nb_wrinkles):
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~WRINKLES~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('WID(s)       : Defect0\n')
		if i < 10:
			parameters.write('Wsize'+str(i)+'(f)    : 1.0\n')  # Amplitude max
			# parameters.write('Wsize'+str(i)+'(f)    : '+str(random.uniform(minWsize, maxWsize))+'\n') # Amplitude max
			parameters.write('Wpos'+str(i)+'(f)     : '+str(random.uniform(minWposX, maxWposX))+','+str(
				random.uniform(minWposY, maxWposY))+','+str(random.uniform(minWposZ, maxWposZ))+'\n')  # center
			# Orientation in degree
			parameters.write('Wori'+str(i)+'(f)     : ' +
			                 str(random.uniform(minWori, maxWori))+'\n')
			parameters.write('Wdamp'+str(i)+'(f)    : '+str(random.uniform(minWdampX, maxWdampX))+','+str(random.uniform(minWdampY,
			                 maxWdampY))+','+str(random.uniform(minWdampZ, maxWdampZ))+'\n')  # reduction of the amplitude through each direction
		else:
			parameters.write('Wsize'+str(i)+'(f)   : ' +
			                 str(random.uniform(minWsize, maxWsize))+'\n')  # Amplitude max
			parameters.write('Wpos'+str(i)+'(f)    : '+str(random.uniform(minWposX, maxWposX))+','+str(
				random.uniform(minWposY, maxWposY))+','+str(random.uniform(minWposZ, maxWposZ))+'\n')  # center
			# Orientation in degree
			parameters.write('Wori'+str(i)+'(f)    : ' +
			                 str(random.uniform(minWori, maxWori))+'\n')
			parameters.write('Wdamp'+str(i)+'(f)   : '+str(random.uniform(minWdampX, maxWdampX))+','+str(random.uniform(minWdampY,
			                 maxWdampY))+','+str(random.uniform(minWdampZ, maxWdampZ))+'\n')  # reduction of the amplitude through each direction
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Abaqus_output(b)  : 1\n') #
	parameters.write('Dune_output(b)    : 0\n') #
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~ABAQUS~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Path2result(s)    : Abaqus/results/\n')
	parameters.write('AbaqusOdbName(s)  : model\n')

	parameters.close()

if __name__ == '__main__':

	write_parameters()
