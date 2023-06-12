import sys
import random
import numpy as np

# geometry in millimeters

def write_parameters():
	# nb_plies = 7
	nb_plies = 2
	nb_wrinkles = 0
	# Xlenght = 500.0
	# Xlenght = 50.0
	Xlenght = 30.0
	# Ylenght = 50.0
	# Ylenght = 50.0
	Ylenght = 5
	# Ylenght = 1.5
	Zlenght = 300.0
	# Zlenght = 100.0
	# Zlenght = 30.0
	height = 0.0

	epsilon = 0.1

	# ply_thickness = Ylenght/(nb_plies+0.0)
	# ply_type = [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0]
	# ply_thickness = [0.15, 0.05, 0.01, 0.05, 0.10, 0.05, 0.01, 0.05, 0.15]
	# ply_thickness = [0.17, 0.03, 0.01, 0.03, 0.14, 0.03, 0.01, 0.03, 0.17]
	# ply_thickness = [0.10, 0.05, 0.025, 0.01, 0.0175, 0.035, 0.07, 0.035, 0.0175, 0.01, 0.025, 0.05, 0.10]
	# tot_ply_thickness = 0.0
	# for t in ply_thickness:
	# 	tot_ply_thickness += t
	# StackSeq = [0.0, 0.0, 0.0, 45.0, 45.0, 45.0, 0.0, -45.0, -45.0]
	# StackSeq = [0.0, 0.0, 0.0, 0.0, 45.0, 45.0, 45.0, 45.0, 45.0, 0.0, -45.0,-45.0, -45.0]

	# StackSeq = [45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 90.0, 0.0, 90.0, 0.0, 90.0, 0.0, -45.0, 45.0, -45.0, 0.0, 45.0, -45.0, 45.0]


	# StackSeq = [0., 0., 0.]
	# StackSeq = [0., 0.]
	
	# large model paper off/on
	# StackSeq = [90.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0]
	# StackSeq = [-45.0, 45.0, 0.0, 0.0, 0.0, 45.0, -45.0]

	# StackSeq = [0., 0., 90., 0.]
	StackSeq = [ 0.0, 0.0 ]
	
	# tot_ply_thickness = Ylenght
	ply_thickness = np.full_like(StackSeq, Ylenght/(nb_plies+0.0))
	# 0 - plies ; 1 - resin // Here: only plies, if resin, it is done via auto resin otherwize need to be specified
	ply_type = np.full_like(StackSeq, 0)

	tot_ply_thickness = Ylenght
	# tot_ply_thickness = Ylenght + epsilon
	# tot_ply_thickness = Ylenght + 2*epsilon
	ply_thickness = np.full_like(StackSeq, Ylenght/(nb_plies+0.0))
	# ply_thickness[0] = epsilon
	ply_type[0] = 2
	# ply_thickness[1] = epsilon
	ply_type[1] = 1

	# ply_thickness[2] = epsilon
	# ply_type[2] = 1
	# ply_thickness[5] = epsilon
	# ply_type[5] = 2


	# WRINKLES Parameters
	minWsize = 1.0
	maxWsize = 1.0
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
	parameters.write('name(s)                  : BeamWeak\n') # mesh name
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
	parameters.write('lc(f)      : 1\n')        # mesh caracteristic size
	parameters.write('dx(i)      : 15\n')       # 15 // 40
	parameters.write('ddy(i)     : 2\n')        #
	# parameters.write('dz(i)      : 160\n')       #
	parameters.write('dz(i)      : 150\n')  # 5,15,15,5 // 10,40,40,10
	parameters.write('dc(i)      : 0\n')        #
	parameters.write('dflange(i) : 0\n')      # discretization of the flange
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~TRANSFORMATION~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('wrinkle(i)        : '+str(nb_wrinkles)+'\n') # Do we need to add wrinkle in the gridMod (c++) code, 2+ for multiple wrinkles
	if nb_wrinkles>0:
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~WRINKLES~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	for i in range(nb_wrinkles):
		parameters.write('WID(s)       : Defect0\n')
		if i<10:
			parameters.write('Wsize'+str(i)+'(f)    : 1.0\n') # Amplitude max
			# parameters.write('Wsize'+str(i)+'(f)    : '+str(random.uniform(minWsize, maxWsize))+'\n') # Amplitude max
			parameters.write('Wpos'+str(i)+'(f)     : '+str(random.uniform(minWposX, maxWposX))+','+str(random.uniform(minWposY, maxWposY))+','+str(random.uniform(minWposZ, maxWposZ))+'\n') # center
			parameters.write('Wori'+str(i)+'(f)     : '+str(random.uniform(minWori, maxWori))+'\n') # Orientation in degree
			parameters.write('Wdamp'+str(i)+'(f)    : '+str(random.uniform(minWdampX, maxWdampX))+','+str(random.uniform(minWdampY, maxWdampY))+','+str(random.uniform(minWdampZ, maxWdampZ))+'\n') # reduction of the amplitude through each direction
		else:
			parameters.write('Wsize'+str(i)+'(f)   : '+str(random.uniform(minWsize, maxWsize))+'\n') # Amplitude max
			parameters.write('Wpos'+str(i)+'(f)    : '+str(random.uniform(minWposX, maxWposX))+','+str(random.uniform(minWposY, maxWposY))+','+str(random.uniform(minWposZ, maxWposZ))+'\n') # center
			parameters.write('Wori'+str(i)+'(f)    : '+str(random.uniform(minWori, maxWori))+'\n') # Orientation in degree
			parameters.write('Wdamp'+str(i)+'(f)   : '+str(random.uniform(minWdampX, maxWdampX))+','+str(random.uniform(minWdampY, maxWdampY))+','+str(random.uniform(minWdampZ, maxWdampZ))+'\n') # reduction of the amplitude through each direction
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~GEO-TRANSFORMATION~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Ramp(b)                : 0\n')
	parameters.write('RotateRVE(b)           : 0\n')
	parameters.write('RotateAxis(b)          : X\n') # "X" or "Z"
	parameters.write('Rotate_start_end(f)    : 100,200\n') # to be matched with dz changes
	parameters.write('AngleRotateRVE(f)      : 90\n') # positive angle
	parameters.write('Rsize(f)               : 6.25\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Abaqus_output(b)  : 0\n')
	parameters.write('Dune_output(b)    : 1\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~ABAQUS~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Path2result(s)         : Abaqus/results/\n')
	parameters.write('AbaqusOdbName(s)       : model\n')
	parameters.write('writeTransformedMSH(b) : 1\n')


	parameters.close()

if __name__ == '__main__':

	write_parameters()
