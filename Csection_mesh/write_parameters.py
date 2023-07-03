import matplotlib.cm as cm
import matplotlib.pyplot as plt
import random
import sys
import numpy as np
np.random.seed(123)

# from skopt.space import Space
# from skopt.sampler import Lhs

# geometry in millimeters

def write_parameters(cnt=-1,p1=0, p2=0, p3=0, p4=0, t_ply=0.196, Zlength = 420.0, height = 55.0, model_name = 'CSpar'):
	# t_ply = ply thickness
	# Zlength = effective length of the spar (between end blocks)
	# total_height = distance from tip of the flanges to the outer surface of the web
	# model_name is the name used for all of the Abaqus input files
	
	nb_plies = 24
	nb_wrinkles = 0
	Xlength = 140.0 # This is the distance along the web of the outer section of the spar between the radii. 140mm is correct for 5mm internal radii
	Ylength = nb_plies*t_ply # laminate thickness
	start = Zlength/2.0 - 150.0 # assumes the end blocks have been placed such that the feature is exactly central
	r_int = 5.0 # internal radius
	height = 55.0 - Ylength - r_int # flange length
	is_coh = 0
	is_GaussianThickness = 0
	is_CornerThickness = 0
	
	# StackSeq = [45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 90.0, 0.0, 90.0, 0.0, 90.0, 0.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0]
	# Jean's model seems to use a left-handed convention for the ply orientations for some reason. Have to swap the -45s and 45s for this to match the right
	# - handed convention used (I think) in the actually specimen
	StackSeq = [-45.0, 45.0, -45.0, 45.0, -45.0, 45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 90.0, 0.0, 90.0, 0.0, 90.0, 0.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0]

	
	ply_thickness = np.full_like(StackSeq, t_ply) # duplicate ply thickness into an array with the same dimensions as StackSeq
	ply_type = np.full_like(StackSeq, 0) # Only plies, resin is done via auto resin
	tot_ply_thickness = Ylength


	# WRINKLES Parameters
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

	if (cnt > 0):
		parameters = open('parameters_'+str(cnt)+'.txt', 'w')
	else:
		parameters = open('parameters.txt', 'w')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('name(s)                  : ' + model_name +'\n') # mesh name
	# parameters.write('name(s)                : '+str(int(Zlenght))+'mm_Cspar\n') # mesh name
	parameters.write('Shape(i)                 : 0\n')   # 0(default): Csection ; 1: flat laminate
	parameters.write('Auto_Resin_betw_plies(b) : 0\n')   # 1: yes ; 0: no (Do I want this to be on?????)
	parameters.write('cohezive_elements(b)     : '+str(is_coh)+'\n')   # 1: yes ; 0: no
	parameters.write('recombine(b)             : 1\n')   # 1: recombine mesh: hex ;  0: no: only prisms
	parameters.write('Shell(b)                 : 0\n')   # 1: recombine mesh: hex ;  0: no: only prisms
	parameters.write('GaussianThickness(b)     : '+str(is_GaussianThickness)+'\n')   # 1: gridtrans ;  0: flat
	parameters.write('CornerThickness(b)       : '+str(is_CornerThickness)+'\n')   # 1: variation ;  0: flat
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	# 6, 12 or 24 # TODO: find a clever way of setting stacking sequence
	parameters.write('np(i)         : '+str(nb_plies)+'\n')
	parameters.write('X(f)          : '+str(Xlength)+'\n')   #
	parameters.write('Y(f)          : '+str(tot_ply_thickness)+'\n')  #
	parameters.write('R(f)          : '+str(r_int)+'\n')    #
	parameters.write('Height(f)     : '+str(height)+'\n')    #
	parameters.write('Z(f)          : '+str(Zlength)+'\n')   #
	parameters.write('e(f)          : 0.01\n')  # Resin layer thickness
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~STACKSEQ~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	# Note - the below if statement is just for neat formatting of the parameter file. Presumably also gmsh defines the coordinate system at 90 degrees to the way we do. Check this
	for i in range(nb_plies):
		if i < 10:
			parameters.write('p'+str(i)+'(f,f,b)   : ' +
			                 str(StackSeq[i]+90.0)+','+str(ply_thickness[i])+','+str(ply_type[i])+'\n')
		else:
			parameters.write('p'+str(i)+'(f,f,b)  : ' +
			                 str(StackSeq[i]+90.0)+','+str(ply_thickness[i])+','+str(ply_type[i])+'\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~MESH~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('lc(f)      : 1\n')    #1 mesh caracteristic size
	parameters.write('dx(i)      : 29\n')   #50 discretization in X direction (plus one)
	parameters.write('dflange(i) : 10\n')    #16 discretization of the flange (plus one)
	parameters.write('ddy(i)     : 2\n')    #2 discretization of ply thickness (plus one)
	parameters.write('dz(i)      : 84\n')   #150 discretization in Z direction (must be a multiple of 84 for ramps to be in the correct place for the 420mm long spar).
	parameters.write('dc(i)      : 5\n')    #12 discretization of corners (plus one)
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~TRANSFORMATION~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	# Do we need to add wrinkle in the gridMod (c++) code, 2+ for multiple wrinkles
	parameters.write('wrinkle(i)        : '+str(nb_wrinkles)+'\n')
	parameters.write('Abaqus_output(b)  : 1\n')
	parameters.write('Dune_output(b)    : 0\n')
	if (nb_wrinkles > 0):
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~WRINKLES~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('WID(s)    : Defect_'+str(cnt)+'\n')
		# parameters.write('Wsize(f)    : '+str(p1)+'\n') # Amplitude max
		# parameters.write('Wpos(f)     : 84.0,73.3,207.7\n') # center
		# parameters.write('Wori(f)     : '+str(p4)+'\n') # Orientation in degree
		# parameters.write('Wdamp(f)    : 8.0,'+str(p2)+','+str(p3)+'\n') # reduction of the amplitude through each direction
		for i in range(nb_wrinkles):
			if i < 10:
				parameters.write('Wsize'+str(i)+'(f)    : ' +
								str(random.uniform(minWsize, maxWsize))+'\n')  # Amplitude max
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
	parameters.write('~~~~~~~~~~~~~~~~GEO-TRANSFORMATION~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Ramp(b)                : 1\n')
	parameters.write('RotateFlanges(b)       : 1\n')
	parameters.write('RotateRVE(b)           : 0\n')
	parameters.write('RotateAxis(b)          : X\n') # "X" or "Z"
	parameters.write('Rotate_start_end(f)    : 100,200\n') # to be matched with dz changes
	parameters.write('AngleRotateRVE(f)      : 90\n') # positive angle
	parameters.write('AngleRotateFlangeRi(f) : 0\n') # positive angle
	parameters.write('AngleRotateFlangeLe(f) : 0\n') # positive angle
	parameters.write('Rsize(f)               : 6.25\n')
	parameters.write('StartEndinZdir(f)      : '+str(start) + ','+str(125.0)+','+str(50.0)+'\n')
	if (is_GaussianThickness):
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~GAUSSIAN~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('Sigma(d)          : 0.1\n')
		parameters.write('Length(d)         : 10\n')
	if (is_CornerThickness):
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~CORNER THICKNESS~~~~~~~~~~~~~~~~~~\n')
		parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		parameters.write('ThicknessVar(d)   : 1\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~ABAQUS~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Path2result(s)    : Abaqus/\n')
	# parameters.write('AbaqusOdbName(s)  : C-spar_FINE_'+str(nb_plies)+'plies\n')
	parameters.write('AbaqusOdbName(s)  : ' + model_name + '\n')

	parameters.close()


def plot_S_dampY(x, title, colors):
    fig, ax = plt.subplots()
    for datax, datay, color in zip(np.array(x)[:, 0], np.array(x)[:, 1], colors):
        plt.scatter(datax, datay, color=color)
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 1], 'bo', label='samples')
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 1], 'bo', markersize=80, alpha=0.5)
    # ax.legend(loc="best", numpoints=1)
    ax.set_xlabel("Size")
    ax.set_xlim([2., 12.])
    ax.set_ylabel("DampY")
    ax.set_ylim([1., 6.])
    plt.title(title)
    plt.savefig(title+".png", dpi=150)


def plot_S_dampZ(x, title, colors):
    fig, ax = plt.subplots()
    for datax, datay, color in zip(np.array(x)[:, 0], np.array(x)[:, 2], colors):
        plt.scatter(datax, datay, color=color)
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 2], 'bo', label='samples')
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 1], 'bo', markersize=80, alpha=0.5)
    # ax.legend(loc="best", numpoints=1)
    ax.set_xlabel("Size")
    ax.set_xlim([2., 12.])
    ax.set_ylabel("DampZ")
    ax.set_ylim([3., 10.])
    plt.title(title)
    plt.savefig(title+".png", dpi=150)


def plot_S_Angle(x, title, colors):
    fig, ax = plt.subplots()
    for datax, datay, color in zip(np.array(x)[:, 0], np.array(x)[:, 3], colors):
        plt.scatter(datax, datay, color=color)
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 3], 'bo', label='samples')
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 1], 'bo', markersize=80, alpha=0.5)
    # ax.legend(loc="best", numpoints=1)
    ax.set_xlabel("Size")
    ax.set_xlim([2., 12.])
    ax.set_ylabel("Angle")
    ax.set_ylim([-5., 5.])
    plt.title(title)
    plt.savefig(title+".png", dpi=150)


def plot_DampYZ(x, title, colors):
    fig, ax = plt.subplots()
    for datax, datay, color in zip(np.array(x)[:, 1], np.array(x)[:, 2], colors):
        plt.scatter(datax, datay, color=color)
    # plt.plot(np.array(x)[:, 1], np.array(x)[:, 2], 'bo', label='samples')
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 1], 'bo', markersize=80, alpha=0.5)
    # ax.legend(loc="best", numpoints=1)
    ax.set_xlabel("DampY")
    ax.set_xlim([1., 6.])
    ax.set_ylabel("DampZ")
    ax.set_ylim([3., 10.])
    plt.title(title)
    plt.savefig(title+".png", dpi=150)


if __name__ == '__main__':

	# EXT 
	a = 40.581
	b = 0.73325

	c = -1.6413
	d = -5.12753

	print((a-b)/(c-d))

	# INT
	a = 34.5916
	b = 0.485142

	c = 5.79578
	d = 2.80836
	print((a-b)/(c-d))

	




	# For Chensen
	# n_samples = 100
	# space = Space([(2., 12.), (1., 6.), (3., 10.), (-5., 5.)])
	# lhs = Lhs(lhs_type="classic", criterion=None)
	# x = lhs.generate(space.dimensions, n_samples)

	# colors = cm.coolwarm(np.linspace(0, 1, n_samples))
	# plot_S_dampY(x, 'Size vs dampY',colors)
	# plot_S_dampZ(x, 'Size vs dampZ',colors)
	# plot_S_Angle(x, 'Size vs Angle',colors)
	# plot_DampYZ(x, 'dampY vs dampZ',colors)

	# f=open('LhsWrinkles.txt', 'w')
	# f.write('Size , DampY, DampZ, Angle\n')
	# for i in range(n_samples):
	# 	f.write(str(np.array(x)[i, 0])+', '+str(np.array(x)[i, 1])+', '+str(np.array(x)[i, 2])+', '+str(np.array(x)[i, 3])+'\n')
	# f.close()

	# for cnt in range(n_samples):
	# 	write_parameters(cnt,np.array(x)[cnt, 0], np.array(x)[cnt, 1], np.array(x)[cnt, 2], np.array(x)[cnt, 3])

	write_parameters()
