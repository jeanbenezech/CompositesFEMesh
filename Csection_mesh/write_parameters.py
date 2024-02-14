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
	Zlength = 740.0  # effective length of the spar (between end blocks)
	# Zlength = 620.0  # effective length of the spar (between end blocks)
	epsilon = 0.01

	# assumes the end blocks have been placed such that the feature is exactly central
	start = Zlength/2.0 - 150.0
	r_int = 5.0  # internal radius
	height = 61.0 - Ylength - r_int  # flange length
	is_coh = 0
	is_GaussianThickness = 0

	## Full stack L2C
	# StackSeq = [ 45.0,  45.0,   0.0,   0.0,
	#      		-45.0, -45.0,  90.0,  90.0, 
	# 			 45.0,  45.0,   0.0,   0.0, 
	# 			-45.0, -45.0,  90.0,  90.0,
	# 			 45.0,  45.0,   0.0,   0.0, 
	# 			-45.0, -45.0,  90.0,  90.0,
	# 			 90.0,  90.0, -45.0, -45.0,
	# 			  0.0,   0.0,  45.0,  45.0,
	# 			 90.0,  90.0, -45.0, -45.0,
	# 			  0.0,   0.0,  45.0,  45.0,
	# 			 90.0,  90.0, -45.0, -45.0,
	# 			  0.0,   0.0,  45.0,  45.0]
	
	###### Ply types: 
		##### 0 : isotropic
		##### 1 : orthotropic composite (ply UD)
		##### 2 : cohesive
		##### 3 : Laminated composite
		##### ## 4 : BC steel (done automaticaly for the C-spar example)

	### no cohesive zone example
	# nb_plies = 24
	# tot_ply_thickness = Ylength
	# ply_t = Ylength/float(nb_plies)
	# ply_thickness = np.full_like(StackSeq, ply_t)
	# ply_type = np.full_like(StackSeq, 1) ## every ply is orthotropic composite UD

	### cohesive zone example
	# nb_plies = 3
	# nb_plies_virtual = 24
	# tot_ply_thickness = Ylength + epsilon
	# ply_t = Ylength/float(nb_plies-1) # nb_plies - #interlayer
	# ply_thickness = np.full_like(StackSeq, ply_t)
	# ply_thickness[2] = epsilon
	# ply_type = np.full_like(StackSeq, 1)
	# ply_type[2] = 2

	### –––––––––– CerTest example 1 ––––––––––
	#### Double plies stack L2C
	nb_plies = 28
	nb_plies_virtual = 24
	########## INNER SURFACE
	StackSeq = [ 45.0, -45.0,
	     		 90.0,  0.0, 0.0, 0.0,
				 45.0, -45.0,
	     		 90.0,   0.0, 
				 45.0, -45.0,
	     		 90.0,   0.0, 
				  0.0,  90.0,
				-45.0,  45.0,
				  0.0,  90.0,
				-45.0,  45.0, 0.0,
				  0.0, 0.0, 90.0,
				-45.0, 45.0]
	########### OUTTER SURFACE
	tot_ply_thickness = Ylength + 4*epsilon
	ply_t = Ylength/float(nb_plies_virtual) # nb_plies - #interlayer
	ply_thickness = np.full_like(StackSeq, ply_t)
	ply_thickness[3] = epsilon
	ply_thickness[5] = epsilon
	ply_thickness[22] = epsilon
	ply_thickness[24] = epsilon
	ply_type = np.full_like(StackSeq, 1)
	ply_type[3] = 2 # CZ
	ply_type[5] = 5 # CZ
	ply_type[22] = 6 # CZ
	ply_type[24] = 7 # CZ


	## –––––––––– CerTest example 2 ––––––––––
	# ## Double plies stack L2C
	# nb_plies = 22
	# nb_plies_virtual = 24
	# ########## INNER SURFACE
	# StackSeq = [ 45.0, -45.0,
	# 	 		 90.0,  0.0, 0.0, 0.0,
	# 			 45.0, -45.0,
	# 	 		 90.0,   0.0, 
	# 			 0.0,
	# 			 0.0,
	# 			  0.0,  90.0,
	# 			-45.0,  45.0, 0.0,
	# 			  0.0, 0.0, 90.0,
	# 			-45.0, 45.0]
	# ########### OUTTER SURFACE
	# tot_ply_thickness = Ylength + 4*epsilon
	# ply_t = Ylength/float(nb_plies_virtual) # nb_plies - #interlayer
	# ply_thickness = np.full_like(StackSeq, ply_t)
	# ply_thickness[3] = epsilon
	# ply_thickness[5] = epsilon
	# ply_thickness[10] = 4*ply_t
	# ply_thickness[11] = 4*ply_t
	# ply_thickness[16] = epsilon
	# ply_thickness[18] = epsilon
	# ply_type = np.full_like(StackSeq, 1)
	# ply_type[3] = 2 # CZ
	# ply_type[5] = 5 # CZ
	# ply_type[10] = 3 # laminate
	# ply_type[11] = 3 # laminate
	# ply_type[16] = 6 # CZ
	# ply_type[18] = 7 # CZ


	### –––––––––– Test case ––––––––––
	# nb_plies = 5
	# nb_plies_virtual = 3
	# # StackSeq = [ 90.0, 0.0, 45.0, 0.0,
	# # 			 -45.0]
	# StackSeq = [ 0.0, 0.0, 0.0, 0.0,
	# 			 0.0]
	# tot_ply_thickness = Ylength + 4*epsilon
	# ply_t = Ylength/float(nb_plies_virtual) # nb_plies - #interlayer
	# ply_thickness = np.full_like(StackSeq, ply_t)
	# ply_thickness[1] = epsilon
	# ply_thickness[3] = epsilon
	# ply_thickness[0] = ply_t*(1.+1./3);
	# ply_thickness[2] = ply_t*(1.+1./3);
	# ply_thickness[4] = ply_t*(1.-2./3);
	# ply_type = np.full_like(StackSeq, 1)
	# ply_type[1] = 2 # CZ
	# ply_type[3] = 5 # CZ

	config = configparser.ConfigParser()
	config.optionxform = lambda option: option ### Preserve capital letters
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GENERAL'] = {}
	config['GENERAL']['name(s)'] = 'CerTest_realistic' ### CerTest, CerTest_realistic, KubernetesTest,CerTest_test,CerTest_test_realistic, CerTest_testflanges, CerTest_testCornerRadius...
	config['GENERAL']['Shape(i)'] = '2' # 1: Laminate type model ; 0: Cspar ; 2 Cspar with specific position of node for BC
	config['GENERAL']['Auto_Resin_betw_plies(b)'] = '0' # add an automatic resin interlayer between plies
	config['GENERAL']['cohezive_elements(b)'] = str(is_coh) # bool have cohesive element
	config['GENERAL']['recombine(b)'] = '1' # 1: recombine mesh: hex ;  0: no: only prisms
	config['GENERAL']['Shell(b)'] = '0' # Modify direction, unique element in thickness
	config['GENERAL']['GaussianThickness(b)'] = '0' # 1: gridtrans ;  0: flat
	config['GENERAL']['CornerThickness(b)'] = '1' # 1: variation ;  0: flat
	config['GENERAL']['Dune_output(b)'] = '1' # write Dune inputs
	config['GENERAL']['Abaqus_output(b)'] = '0' # write Abaqus inputs
	config['GENERAL']['Path2result(s)'] = 'Abaqus/results/' # write Abaqus inputs
	config['GENERAL']['AbaqusOdbName(s)'] = 'CSPAR' # write Abaqus inputs
	config['GENERAL']['writeTransformedMSH(b)'] = '0' # write Abaqus grid with transformation

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	config['GEOMETRY'] = {}
	config['GEOMETRY']['np(i)'] = str(nb_plies) # number of plies
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
	config['MESH']['dx(i)'] = '110' # 25 ;; 110 # Discretization top part of the C
	config['MESH']['dflange(i)'] = '55' # 55 # uniform discretization of the flange for shape = 1
	config['MESH']['ddy(i)'] = '2' # 2 # ply discretization (number of DoF used so 2 = 1 element)
	config['MESH']['dz(i)'] = '1,20,90,20,1' # 150 // 1,9,25,9,1 ;; 1,20,90,20,1 # Discretization of the part lenght. if 1 int is provided: uniform; else multiple discretization
	config['MESH']['dc(i)'] = '15' # 7 ;; 15 # Discretization the corner
	config['MESH']['dflange_bot(i)'] = '12' # 6 ;; 12 # Discretization the flange bellow BC node
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
	# config['GEO-TRANSFORMATION']['dZInterval(f)'] = '100,260,360,520' # Start-stop of each section
	config['GEO-TRANSFORMATION']['dZInterval(f)'] = '160,320,420,580' # Start-stop of each section
	config['GEO-TRANSFORMATION']['dFlangeInterval(b)'] = '36.27,38.27,40.27,42.27,44.27,46.27' # How many sections with different discretization we have for flanges. Each layer created will be 1 element thick
	config['GEO-TRANSFORMATION']['RotateFlanges(b)'] = '1' # Allows flanges rotation
	config['GEO-TRANSFORMATION']['RotateRVE(b)'] = '0' # Rotation of the RVE in the RVE case
	config['GEO-TRANSFORMATION']['Rotate_start_end(f)'] = '100,200' # Strat-stop rotation to create corner test case
	config['GEO-TRANSFORMATION']['AngleRotateRVE(f)'] = '90' # positive angle rotation
	config['GEO-TRANSFORMATION']['RotateAxis(s)'] = 'X' # "X" or "Z"
	config['GEO-TRANSFORMATION']['AngleRotateFlangeRi(f)'] = '1.0' # positive angle
	config['GEO-TRANSFORMATION']['AngleRotateFlangeLe(f)'] = '0.5' # positive angle
	config['GEO-TRANSFORMATION']['ThicknessVar(d)'] = '1.' # ratio thickness of the corners

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
