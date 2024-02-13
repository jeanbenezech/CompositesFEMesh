import numpy as np

import sys
import random
import matplotlib.pyplot as plt
# from skopt.space import Space
# from skopt.sampler import Lhs
import matplotlib.cm as cm
import configparser

# geometry in millimeters

def write_parameters(cnt=-1,p1=0, p2=0, p3=0, p4=0):
	nb_plies = 1 #fixed
	nb_wrinkles = 0 # let it to 0
	Xlength = 140.0
	Ylength = 5.0 #6.48
	Zlength = 420.0
	height = 35 # height of the flanges

	# assumes the end blocks have been placed such that the feature is exactly central
	start = Zlength/2.0 - 150.0
	r_int = 5.0  # internal radius
	height = 50.0 - Ylength - r_int  # flange length
	r_int = 5.0 # internal radius
	is_coh = 0 # Not for shell

	is_GaussianThickness = 1
	is_CornerThickness = 1

	# ply_thickness = Ylenght/(nb_plies+0.0)
	StackSeq = [0.0]

	ply_thickness = np.full_like(StackSeq, Ylength/(nb_plies+0.0))
	ply_type = np.full_like(StackSeq, 0) # Only plies, resin is done via auto resin


	config = configparser.ConfigParser()
	config.optionxform = lambda option: option ### Preserve capital letters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GENERAL'] = {}
	config['GENERAL']['name(s)'] = 'C-spar_continuum_shell' ### CerTest,KubernetesTest,CerTest_test...
	config['GENERAL']['Shape(i)'] = '0' # 1: Laminate type model ; 0: Cspar ; 2 Cspar with specific position of node for BC
	config['GENERAL']['Auto_Resin_betw_plies(b)'] = '0' # add an automatic resin interlayer between plies
	config['GENERAL']['cohezive_elements(b)'] = str(is_coh) # bool have cohesive element
	config['GENERAL']['recombine(b)'] = '1' # 1: recombine mesh: hex ;  0: no: only prisms
	config['GENERAL']['Shell(b)'] = '1' # Modify direction, unique element in thickness
	config['GENERAL']['GaussianThickness(b)'] = str(is_GaussianThickness) # 1: gridtrans ;  0: flat
	config['GENERAL']['CornerThickness(b)'] = str(is_CornerThickness) # 1: variation ;  0: flat
	config['GENERAL']['Dune_output(b)'] = '0' # write Dune inputs
	config['GENERAL']['Abaqus_output(b)'] = '1' # write Abaqus inputs
	config['GENERAL']['Path2result(s)'] = 'Abaqus/results/' # write Abaqus inputs
	config['GENERAL']['AbaqusOdbName(s)'] = 'C-spar_continuum_shell' # write Abaqus inputs
	config['GENERAL']['writeTransformedMSH(b)'] = '0' # write Abaqus grid with transformation

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GEOMETRY'] = {}
	config['GEOMETRY']['np(i)'] = str(nb_plies) # number of plies
	config['GEOMETRY']['X(f)'] = str(Xlength) # Beam: witdh; Cspar: width of the ...
	config['GEOMETRY']['Y(f)'] = str(Ylength) # number of plies
	config['GEOMETRY']['R(f)'] = str(r_int) # internal radius of the corners of the C shape; not used in laminate
	config['GEOMETRY']['Height(f)'] = str(height) # height of the flanges from bottom to the corner
	config['GEOMETRY']['Z(f)'] = str(Zlength) # Lenght of the spar of beam
	config['GEOMETRY']['e(f)'] = '0.01' # interlayer thickness

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~STACKSEQ~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['STACKSEQ'] = {}
	for i in range(nb_plies):
		config['STACKSEQ']['p'+str(i)+'(f,f,b)'] = str(StackSeq[i]+90.0)+','+str(ply_thickness[i])+','+str(int(ply_type[i]))

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~MESH~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['MESH'] = {}
	config['MESH']['lc(f)'] = '1' # 1 mesh caracteristic size, overwritten if any following parameter is used
	config['MESH']['dx(i)'] = '10' # 100 # Discretization top part of the C
	config['MESH']['dflange(i)'] = '4' # 55 # uniform discretization of the flange for shape = 1
	config['MESH']['ddy(i)'] = '2' # 2 # ply discretization (number of DoF used so 2 = 1 element)
	config['MESH']['dz(i)'] = '40' # 150 // 1,9,25,9,1 # Discretization of the part lenght. if 1 int is provided: uniform; else multiple discretization
	config['MESH']['dc(i)'] = '6' # 15 # Discretization the corner
	config['MESH']['dflange_bot(i)'] = '6' # 12 # Discretization the flange bellow BC node
	config['MESH']['dflange_top(i)'] = '2' # 12 # Discretization the flange above BC node

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~GEO-TRANSFORMATION~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GEO-TRANSFORMATION'] = {}
	config['GEO-TRANSFORMATION']['Ramp(b)'] = '1' # Add a ramp in the C-spar example
	config['GEO-TRANSFORMATION']['Rsize(f)'] = '6.25' # Ramp size
	config['GEO-TRANSFORMATION']['StartEndinZdir(f)'] = str(start) +','+str(125.0)+','+str(50.0) # Starts-stops of the ramp
	config['GEO-TRANSFORMATION']['RigidBoundary(b)'] = '0' # Layer of a rigid boundary for simply supported BC
	config['GEO-TRANSFORMATION']['dZIntervalsize(i)'] = '4' # How many sections with different discretization we have in length direction (linked to dz in MESH section)
	config['GEO-TRANSFORMATION']['dZInterval(f)'] = '100,260,360,520' # Start-stop of each section
	config['GEO-TRANSFORMATION']['dFlangeInterval(b)'] = '36.27,38.27,40.27,42.27,44.27' # How many sections with different discretization we have for flanges. Each layer created will be 1 element thick
	config['GEO-TRANSFORMATION']['RotateFlanges(b)'] = '1' # Allows flanges rotation
	config['GEO-TRANSFORMATION']['RotateRVE(b)'] = '0' # Rotation of the RVE in the RVE case
	config['GEO-TRANSFORMATION']['Rotate_start_end(f)'] = '100,200' # Strat-stop rotation to create corner test case
	config['GEO-TRANSFORMATION']['AngleRotateRVE(f)'] = '90' # positive angle rotation
	config['GEO-TRANSFORMATION']['RotateAxis(s)'] = 'X' # "X" or "Z"
	config['GEO-TRANSFORMATION']['AngleRotateFlangeRi(f)'] = '3' # positive angle
	config['GEO-TRANSFORMATION']['AngleRotateFlangeLe(f)'] = '10' # positive angle
	config['GEO-TRANSFORMATION']['ThicknessVar(d)'] = '1' # ratio thickness of the corners

	if (is_GaussianThickness):
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# ~~~~~~~~~~~~~~~~~~~~GAUSSIAN~~~~~~~~~~~~~~~~~~~~~~
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		config['GAUSSIAN'] = {}
		config['GAUSSIAN']['Sigma(d)'] = '0.1' # covariance
		config['GAUSSIAN']['Length(d)'] = '10.0' # mean

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~WRINKLES~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['WRINKLES'] = {}
	config['WRINKLES']['wrinkle(i)'] = str(nb_wrinkles) # bool do we add wrinkle in the mesh
	if (nb_wrinkles > 0):
		### WRINKLES Parameters
		minWsize = -0.2
		maxWsize = 0.2
		minWposX = 0.11*Xlength
		maxWposX = 0.89*Xlength
		minWposY = height - 0.31*Ylength
		maxWposY = height + 0.0
		minWposZ = 0.11*Zlength
		maxWposZ = 0.89*Zlength
		minWdampX = 2.0
		maxWdampX = 9.0
		minWdampY = 0.6
		maxWdampY = 1.0
		minWdampZ = 1.5
		maxWdampZ = 2.5
		minWori = -20
		maxWori = 20
		config['WRINKLES']['WID(s)'] = 'Defect_'+str(cnt) # wrinkle ID
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
