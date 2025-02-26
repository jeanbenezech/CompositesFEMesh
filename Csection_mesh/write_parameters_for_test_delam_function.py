# import matplotlib.cm as cm
# import matplotlib.pyplot as plt
import random
import sys
import numpy as np
# np.random.seed(123)
import configparser

def write_parameters(cnt=-1, p1=0, p2=0, p3=0, p4=0):
	nb_wrinkles = 0
	# geometry in millimeters
	Xlength = 140.0  # This is the distance along the web of the outer section of the spar between the radii. 140mm is correct for 5mm internal radii
	# Ylength = nb_plies*t_ply  # laminate thickness
	Ylength = 6.  # laminate thickness
	# Zlength = 420.0  # effective length of the spar (between end blocks)
	Zlength = 740.0  # effective length of the spar (with 160mm end blocks)
	# Zlength = 420.0  # effective length of the spar (between end blocks)
	# Zlength = 620.0  # effective length of the spar (with 100mm end blocks)
	epsilon = 0.01

	# assumes the end blocks have been placed such that the feature is exactly central
	start = Zlength/2.0 - 150.0
	r_int = 5.0  # internal radius
	height = 61.0 - Ylength - r_int  # flange length
	is_coh = 0
	is_GaussianThickness = 0

	# ### Test model for CZ pre damage zone shape validations
	Interlayer_position = 20 ## 20 will be between layer 20 and 21
	nb_plies = 3
	nb_plies_virtual = 3
	tot_ply_thickness = Ylength + epsilon
	ply_t = Ylength/float(nb_plies_virtual) # nb_plies - #interlayer
	StackSeq = [0.0, 0.0, 0.0]
	ply_thickness = np.full_like(StackSeq, ply_t)
	ply_thickness[0] = (Interlayer_position+0.) * Ylength/24.
	ply_thickness[1] = epsilon
	ply_thickness[2] = (24-Interlayer_position+0.) * Ylength/24.
	ply_type = np.full_like(StackSeq, 1)
	ply_type[1] = 2



	config = configparser.ConfigParser()
	config.optionxform = lambda option: option ### Preserve capital letters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GENERAL'] = {}
	config['GENERAL']['name(s)'] = 'Cspar_test_delam' ### CerTest, CerTest_realistic, KubernetesTest,CerTest_test,CerTest_test_realistic, CerTest_testflanges, CerTest_testCornerRadius...
	config['GENERAL']['Shape(i)'] = '2' # 1: Laminate type model ; 0: Cspar ; 2 Cspar with specific position of node for BC
	config['GENERAL']['Auto_Resin_betw_plies(b)'] = '0' # add an automatic resin interlayer between plies
	config['GENERAL']['cohezive_elements(b)'] = str(is_coh) # bool have cohesive element
	config['GENERAL']['Abaqus_cohezive_ids_size(i)'] = '2' # bool have cohesive element
	config['GENERAL']['Abaqus_cohezive_ids(i)'] = '2,4' # bool have cohesive element // in terms of ply number
	config['GENERAL']['recombine(b)'] = '1' # 1: recombine mesh: hex ;  0: no: only prisms
	config['GENERAL']['Shell(b)'] = '0' # Modify direction, unique element in thickness
	config['GENERAL']['GaussianThickness(b)'] = '0' # 1: gridtrans ;  0: flat
	config['GENERAL']['CornerThickness(b)'] = '0' # 1: variation ;  0: flat
	config['GENERAL']['Dune_output(b)'] = '0' # write Dune inputs
	config['GENERAL']['Abaqus_output(b)'] = '0' # write Abaqus inputs
	config['GENERAL']['Path2result(s)'] = 'Abaqus/results/' # write Abaqus inputs
	config['GENERAL']['AbaqusOdbName(s)'] = 'Cohesive_test' # write Abaqus inputs
	config['GENERAL']['writeTransformedMSH(b)'] = '0' # write Abaqus grid with transformation


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~EXTRA PARAM~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['TESTDELAM'] = {}
	config['TESTDELAM']['test_delam(b)'] = '1' # [0 or 1] write Abaqus grid with transformation
	config['TESTDELAM']['flatten(b)'] = '1' # [0 or 1] write Abaqus grid with transformation
	config['TESTDELAM']['sizeCentersZ(i)'] = '2' # size of vector line below
	config['TESTDELAM']['centersZ(f)'] = '320., 380.' # position of center of delam in Z dir
	config['TESTDELAM']['sizeCentersX(i)'] = '5' # size of vector line below
	config['TESTDELAM']['centersX(f)'] = '-50., 0., 50., 100., 150.' # position of center of delam in curviliear XY dir


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GEOMETRY'] = {}
	config['GEOMETRY']['np(i)'] = str(nb_plies) # number of plies
	config['GEOMETRY']['npreal(i)'] = str(nb_plies_virtual) # number of plies
	config['GEOMETRY']['X(f)'] = str(Xlength) # Beam: witdh; Cspar: width of the ...
	config['GEOMETRY']['Y(f)'] = str(tot_ply_thickness) # number of plies
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
	config['MESH']['dx(i)'] = '110' # 25 ;; 60 ;; 110 # Discretization top part of the C
	config['MESH']['dflange(i)'] = '55' # 55 # uniform discretization of the flange for shape = 1
	config['MESH']['ddy(i)'] = '2' # 2 # ply discretization (number of DoF used so 2 = 1 element)
	######## we discretize in Z dir properly only the central part of the cspar
	config['MESH']['dz(i)'] = '1,1,90,1,1' # 150
	config['MESH']['dc(i)'] = '15' # 7 ;; 15 # Discretization the corner
	config['MESH']['dflange_bot(i)'] = '17' # 6 ;; 12 ;; 17 # Discretization the flange bellow BC node
	config['MESH']['dflange_top(i)'] = '4' # 2 ;; 4 # Discretization the flange above BC node

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~GEO-TRANSFORMATION~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GEO-TRANSFORMATION'] = {}
	config['GEO-TRANSFORMATION']['Ramp(b)'] = '1' # Add a ramp in the C-spar example
	config['GEO-TRANSFORMATION']['Rsize(f)'] = '6.25' # Ramp size
	config['GEO-TRANSFORMATION']['StartEndinZdir(f)'] = str(start) +','+str(125.0)+','+str(50.0) # Starts-stops of the ramp
	config['GEO-TRANSFORMATION']['RigidBoundary(b)'] = '1' # Layer of a rigid boundary for simply supported BC
	config['GEO-TRANSFORMATION']['dZIntervalsize(i)'] = '4' # How many sections with different discretization we have in length direction (linked to dz in MESH section)
	config['GEO-TRANSFORMATION']['dZInterval(f)'] = '160,280,460,580' # Start-stop of each section
	config['GEO-TRANSFORMATION']['dFlangeInterval(b)'] = '38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0' # How many sections with different discretization we have for flanges. Each layer created will be 1 element thick
	config['GEO-TRANSFORMATION']['RotateFlanges(b)'] = '0' # Allows flanges rotation
	config['GEO-TRANSFORMATION']['RotateRVE(b)'] = '0' # Rotation of the RVE in the RVE case
	config['GEO-TRANSFORMATION']['Rotate_start_end(f)'] = '100,200' # Strat-stop rotation to create corner test case
	config['GEO-TRANSFORMATION']['AngleRotateRVE(f)'] = '90' # positive angle rotation
	config['GEO-TRANSFORMATION']['RotateAxis(s)'] = 'X' # "X" or "Z"
	config['GEO-TRANSFORMATION']['AngleRotateFlangeRi(f)'] = '1.0' # positive angle
	config['GEO-TRANSFORMATION']['AngleRotateFlangeLe(f)'] = '0.5' # positive angle
	config['GEO-TRANSFORMATION']['ThicknessVar(d)'] = '0.2' # ratio thickness of the corners

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

	if (cnt > 0):
		with open('parameters_'+str(cnt)+'.ini', 'w') as configfile:
			config.write(configfile)
	else:
		with open('parameters.ini', 'w') as configfile:
			config.write(configfile)


if __name__ == '__main__':

	write_parameters()
