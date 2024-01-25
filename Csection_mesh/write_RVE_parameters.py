import sys
import random
import numpy as np
import configparser

# geometry in millimeters

def write_parameters():
	nb_plies = 24
	nb_wrinkles = 0
	Xlength = 100.0
	Ylength = 50.0 # thickness
	Zlength = 100.0
	Radius = 10.0
	height = 0.0

	StackSeq = [45.0, -45.0, 45.0, -45.0, 
			 	45.0, -45.0, 0.0, 90.0, 
				0.0, 90.0, 0.0, 90.0, 
				90.0, 0.0, 90.0, 0.0, 
				90.0, 0.0, -45.0, 45.0, 
				-45.0, 45.0, -45.0, 45.0]
	tot_ply_thickness = Ylength
	ply_thickness = np.full_like(StackSeq, Ylength/(nb_plies+0.0))
	ply_type = np.full_like(StackSeq, 0) # 0 - plies ; 1 - resin // Here: only plies, if resin, it is done via auto resin otherwize need to be specified



	config = configparser.ConfigParser()
	config.optionxform = lambda option: option ### Preserve capital letters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GENERAL'] = {}
	config['GENERAL']['name(s)'] = 'RVE' ### CerTest,KubernetesTest,CerTest_test...
	config['GENERAL']['Shape(i)'] = '1' # 1: Laminate type model ; 0: Cspar ; 2 Cspar with specific position of node for BC
	config['GENERAL']['Auto_Resin_betw_plies(b)'] = '0' # add an automatic resin interlayer between plies
	config['GENERAL']['cohezive_elements(b)'] = '0' # bool have cohesive element
	config['GENERAL']['recombine(b)'] = '1' # 1: recombine mesh: hex ;  0: no: only prisms
	config['GENERAL']['Shell(b)'] = '0' # Modify direction, unique element in thickness
	config['GENERAL']['GaussianThickness(b)'] = '0' # 1: gridtrans ;  0: flat
	config['GENERAL']['CornerThickness(b)'] = '0' # 1: variation ;  0: flat
	config['GENERAL']['Dune_output(b)'] = '0' # write Dune inputs
	config['GENERAL']['Abaqus_output(b)'] = '1' # write Abaqus inputs
	config['GENERAL']['Path2result(s)'] = 'Abaqus/results/' # write Abaqus inputs
	config['GENERAL']['AbaqusOdbName(s)'] = 'RVE' # write Abaqus inputs
	config['GENERAL']['writeTransformedMSH(b)'] = '1' # write Abaqus grid with transformation

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GEOMETRY'] = {}
	config['GEOMETRY']['np(i)'] = str(nb_plies) # number of plies
	config['GEOMETRY']['X(f)'] = str(Xlength) # Beam: witdh; Cspar: width of the ...
	config['GEOMETRY']['Y(f)'] = str(tot_ply_thickness) # number of plies
	config['GEOMETRY']['Z(f)'] = str(Zlength) # Lenght of the spar of beam
	config['GEOMETRY']['R(f)'] = str(Radius) # Lenght of the spar of beam
	config['GEOMETRY']['e(f)'] = '0.01' # interlayer thickness

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~STACKSEQ~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['STACKSEQ'] = {}
	for i in range(nb_plies):
		config['STACKSEQ']['p'+str(i)+'(f,f,b)'] = str(StackSeq[i])+','+str(ply_thickness[i])+','+str(int(ply_type[i]))


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~MESH~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['MESH'] = {}
	config['MESH']['lc(f)'] = '1' # 1 mesh caracteristic size, overwritten if any following parameter is used
	config['MESH']['dx(i)'] = '20' # 15 # Discretization top part of the C
	config['MESH']['ddy(i)'] = '3' # 2 # ply discretization (number of DoF used so 2 = 1 element)
	config['MESH']['dz(i)'] = '20' # 150 // 5,15,15,5 # Discretization of the part lenght. if 1 int is provided: uniform; else multiple discretization

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~GEO-TRANSFORMATION~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GEO-TRANSFORMATION'] = {}
	config['GEO-TRANSFORMATION']['Ramp(b)'] = '0' # Add a ramp in the C-spar example
	config['GEO-TRANSFORMATION']['Rsize(f)'] = '0.0' # Ramp size
	# config['GEO-TRANSFORMATION']['StartEndinZdir(f)'] = str(start) +','+str(125.0)+','+str(50.0) # Starts-stops of the ramp
	config['GEO-TRANSFORMATION']['RigidBoundary(b)'] = '0' # Layer of a rigid boundary for simply supported BC
	# config['GEO-TRANSFORMATION']['dZIntervalsize(i)'] = '2' # How many sections with different discretization we have in length direction (linked to dz in MESH section)
	# config['GEO-TRANSFORMATION']['dZInterval(f)'] = '5,10' # Start-stop of each section
	config['GEO-TRANSFORMATION']['RotateRVE(b)'] = '1' # Rotation of the RVE in the RVE case
	config['GEO-TRANSFORMATION']['Rotate_start_end(f)'] = '0,100' # Strat-stop rotation to create corner test case
	config['GEO-TRANSFORMATION']['AngleRotateRVE(f)'] = '45.' # positive angle rotation
	# config['GEO-TRANSFORMATION']['InteriorRadiusRVE(f)'] = '100.' # computed from anlgeRVE and Rotate_start_end
	config['GEO-TRANSFORMATION']['RotateAxis(s)'] = 'X' # "X" or "Z"
	# config['GEO-TRANSFORMATION']['dFlangeInterval(b)'] = '36.27,38.27,40.27,42.27,44.27' # How many sections with different discretization we have for flanges. Each layer created will be 1 element thick
	# config['GEO-TRANSFORMATION']['RotateFlanges(b)'] = '0' # Allows flanges rotation
	# config['GEO-TRANSFORMATION']['AngleRotateFlangeRi(f)'] = '20' # positive angle
	# config['GEO-TRANSFORMATION']['AngleRotateFlangeLe(f)'] = '20' # positive angle

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~WRINKLES~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['WRINKLES'] = {}
	config['WRINKLES']['wrinkle(i)'] = str(nb_wrinkles) # bool do we add wrinkle in the mesh
	if (nb_wrinkles > 0):
		### WRINKLES Parameters
		minWsize = 1.0
		maxWsize = 1.0
		minWposX = 0.30*Xlength
		maxWposX = 0.70*Xlength
		minWposY = 0.49* Ylength
		maxWposY = 0.51* Ylength
		minWposZ = 0.45*Zlength
		maxWposZ = 0.55*Zlength
		minWdampX = 0.01
		maxWdampX = 0.06
		minWdampY = 0.05
		maxWdampY = 0.15
		minWdampZ = 0.05
		maxWdampZ = 0.10
		minWori = -20
		maxWori = 20
		config['WRINKLES']['WID(s)'] = 'Defect_0' # wrinkle ID
		for i in range(nb_wrinkles):
			config['WRINKLES']['Wsize'+str(i)+'(f)'] = str(random.uniform(minWsize, maxWsize)) # Amplitude max
			config['WRINKLES']['Wpos'+str(i)+'(f)'] = str(random.uniform(minWposX, maxWposX))+','+\
													  str(random.uniform(minWposY, maxWposY))+','+\
													  str(random.uniform(minWposZ, maxWposZ)) # Center
			config['WRINKLES']['Wori'+str(i)+'(f)'] = str(random.uniform(minWori, maxWori)) # Orientation in degree
			config['WRINKLES']['Wdamp'+str(i)+'(f)'] = str(random.uniform(minWdampX, maxWdampX))+','+\
													  str(random.uniform(minWdampY, maxWdampY))+','+\
													  str(random.uniform(minWdampZ, maxWdampZ)) # reduction of the amplitude through each direction

	with open('parameters.ini', 'w') as configfile:
		config.write(configfile)

if __name__ == '__main__':

	write_parameters()
