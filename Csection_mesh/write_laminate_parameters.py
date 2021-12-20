import sys
import random

# geometry in millimeters

def write_parameters():
	nb_plies = 4
	nb_wrinkles = 0
	Xlenght = 140.0
	Ylenght = 5.00
	Zlenght = 500.0
	height = 0.0

	ply_thickness = Ylenght/(nb_plies+0.0)
	StackSeq = [0.0, 90.0, 0.0, 90.0]

	# WRINKLES Parameters
	minWsize = -0.5
	maxWsize = 0.8

	minWposX = 0.11*Xlenght
	maxWposX = 0.89*Xlenght
	# minWposX = 0.49*Xlenght
	# maxWposX = 0.51*Xlenght

	minWposY = -0.3* Ylenght
	maxWposY = 0#Ylenght

	minWposZ = 0.45*Zlenght
	maxWposZ = 0.55*Zlenght

	minWdampX = 2.0
	maxWdampX = 9.0

	minWdampY = 0.6
	maxWdampY = 1.0

	minWdampZ = 1.5
	maxWdampZ = 2.5

	minWori = -20
	maxWori = 20

	parameters=open('parameters.txt','w')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('name(s)                : beam\n') # mesh name
	parameters.write('Shape(i)               : 1\n')   # 0(default): Csection ; 1: flat laminate
	parameters.write('Resin_betw_plies(b)    : 1\n')   # 1: yes ; 0: no
	parameters.write('cohezive_elements(b)   : 0\n')   # 1: yes ; 0: no
	parameters.write('recombine(b)           : 1\n')   # 1: recombine mesh: hex ;  0: no: only prisms
	parameters.write('Shell(b)               : 0\n')   # 1: recombine mesh: hex ;  0: no: only prisms
	parameters.write('GaussianThickness(b)   : 0\n')   # 1: gridtrans ;  0: flat
	parameters.write('CornerThickness(b)     : 0\n')   # 1: variation ;  0: flat
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
			parameters.write('p'+str(i)+'(f,f)   : '+str(StackSeq[i])+','+str(ply_thickness)+'\n')
		else:
			parameters.write('p'+str(i)+'(f,f)  : '+str(StackSeq[i])+','+str(ply_thickness)+'\n')
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
	parameters.write('dx(i)    : 40\n')       #
	parameters.write('ddy(i)   : 3\n')        #
	parameters.write('dz(i)    : 100\n')       #
	parameters.write('dc(i)    : 0\n')        #
	parameters.write('dflange(i) : 0\n')     # discretization of the flange
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
		parameters.write('WID(s)             : Defect0\n')
		if i<10:
			parameters.write('Wsize'+str(i)+'(f)    : '+str(random.uniform(minWsize, maxWsize))+'\n') # Amplitude max
			parameters.write('Wpos'+str(i)+'(f)     : '+str(random.uniform(minWposX, maxWposX))+','+str(random.uniform(minWposY, maxWposY))+','+str(random.uniform(minWposZ, maxWposZ))+'\n') # center
			parameters.write('Wori'+str(i)+'(f)     : '+str(random.uniform(minWori, maxWori))+'\n') # Orientation in degree
			parameters.write('Wdamp'+str(i)+'(f)    : '+str(random.uniform(minWdampX, maxWdampX))+','+str(random.uniform(minWdampY, maxWdampY))+','+str(random.uniform(minWdampZ, maxWdampZ))+'\n') # reduction of the amplitude through each direction
		else:
			parameters.write('Wsize'+str(i)+'(f)   : '+str(random.uniform(minWsize, maxWsize))+'\n') # Amplitude max
			parameters.write('Wpos'+str(i)+'(f)    : '+str(random.uniform(minWposX, maxWposX))+','+str(random.uniform(minWposY, maxWposY))+','+str(random.uniform(minWposZ, maxWposZ))+'\n') # center
			parameters.write('Wori'+str(i)+'(f)    : '+str(random.uniform(minWori, maxWori))+'\n') # Orientation in degree
			parameters.write('Wdamp'+str(i)+'(f)   : '+str(random.uniform(minWdampX, maxWdampX))+','+str(random.uniform(minWdampY, maxWdampY))+','+str(random.uniform(minWdampZ, maxWdampZ))+'\n') # reduction of the amplitude through each direction
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~RAMP~~~~~~~~~~~~~~~~~~~~~~~~\n')
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
