import sys
import random
import numpy as np

# geometry in millimeters

def write_parameters():
	nb_plies = 2
	nb_wrinkles = 0
	# Xlenght = 500.0
	Xlenght = 200.0
	# Xlenght = 140.0
	# Ylenght = 50.0
	Ylenght = 50.0
	# Zlenght = 600.0
	# Zlenght = 10.0
	Zlenght = 400.0
	height = 0.0

	StackSeq = [0.0, 0.0]
	# StackSeq = [45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 90.0, 0.0, 90.0, 0.0, 90.0, 0.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0]

	ply_thickness = np.full_like(StackSeq, Ylenght/(nb_plies+0.0))
	ply_type = np.full_like(StackSeq, 1) # Only plies, resin is done via auto resin


	# WRINKLES Parameters
	minWsize = 1.0
	maxWsize = 1.0

	minWposX = 0.30*Xlenght
	maxWposX = 0.70*Xlenght
	# minWposX = 0.49*Xlenght
	# maxWposX = 0.51*Xlenght

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
	parameters.write('name(s)                  : nlbeam\n') # mesh name
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
	parameters.write('Y(f)          : '+str(Ylenght)+'\n')  #
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
	# for i in range(8):
	# 	if i<10:
	# 		parameters.write('p'+str(i)+'(f,f)   : 0.0,'+str(ply_thickness)+'\n')
	# 	else:
	# 		parameters.write('p'+str(i)+'(f,f)  : 0.0,'+str(ply_thickness)+'\n')
	# parameters.write('p8(f,f)   : 0.0,0.31\n')
	# for i in range(9,17):
	# 	if i<10:
	# 		parameters.write('p'+str(i)+'(f,f)   : 0.0,'+str(ply_thickness)+'\n')
	# 	else:
	# 		parameters.write('p'+str(i)+'(f,f)  : 0.0,'+str(ply_thickness)+'\n')
	# parameters.write('p17(f,f)  : 0.0,0.31\n')
	# for i in range(18,20):
	# 	if i<10:
	# 		parameters.write('p'+str(i)+'(f,f)   : 0.0,'+str(ply_thickness)+'\n')
	# 	else:
	# 		parameters.write('p'+str(i)+'(f,f)  : 0.0,'+str(ply_thickness)+'\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~MESH~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('lc(f)    : 1\n')        # mesh caracteristic size
	parameters.write('dx(i)    : 7\n')       #
	parameters.write('ddy(i)   : 2\n')        #
	parameters.write('dz(i)    : 15\n')      #
	parameters.write('dc(i)    : 0\n')        #
	parameters.write('dflange(i) : 0\n')      # discretization of the flange
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~TRANSFORMATION~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('wrinkle(i)        : '+str(nb_wrinkles)+'\n') # Do we need to add wrinkle in the gridMod (c++) code, 2+ for multiple wrinkles
	parameters.write('Ramp(b)           : 0\n') #
	parameters.write('Abaqus_output(b)  : 0\n') #
	parameters.write('Dune_output(b)    : 1\n') #
	for i in range(nb_wrinkles):
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~WRINKLES~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
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
	parameters.write('Rsize(f)          : 6.25\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~ABAQUS~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Path2result(s)    : Abaqus/results/\n')
	parameters.write('AbaqusOdbName(s)  : model\n')


	parameters.close()

if __name__ == '__main__':

	write_parameters()
