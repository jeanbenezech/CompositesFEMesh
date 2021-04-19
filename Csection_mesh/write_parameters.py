import numpy as np
np.random.seed(123)

import sys
import random
import matplotlib.pyplot as plt
from skopt.space import Space
from skopt.sampler import Lhs
import matplotlib.cm as cm

# geometry in millimeters

def write_parameters(cnt,p1, p2, p3, p4):
	nb_plies = 6
	nb_wrinkles = 1
	Xlenght = 150.0
	Ylenght = 7.0
	Zlenght = 500.0
	height = 75

	ply_thickness = Ylenght/(nb_plies+0.0)
	StackSeq = [45.0, -45.0, 90.0, 0.0, -45.0, 45.0]

	# WRINKLES Parameters
	minWsize = -0.2
	maxWsize = 0.2

	minWposX = 0.11*Xlenght
	maxWposX = 0.89*Xlenght

	minWposY = height - 0.31*Ylenght
	maxWposY = height + 0.0

	minWposZ = 0.11*Zlenght
	maxWposZ = 0.89*Zlenght

	minWdampX = 2.0
	maxWdampX = 9.0

	minWdampY = 0.6
	maxWdampY = 1.0

	minWdampZ = 1.5
	maxWdampZ = 2.5

	minWori = -20
	maxWori = 20



	parameters=open('parameters_'+str(cnt)+'.txt','w')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('name(s)                : Csection\n') # mesh name
	parameters.write('Shape(i)               : 0\n')   # 0(default): Csection ; 1: flat laminate
	parameters.write('Resin_betw_plies(b)    : 0\n')   # 1: yes ; 0: no
	parameters.write('cohezive_elements(b)   : 0\n')   # 1: yes ; 0: no
	parameters.write('recombine(b)           : 1\n')   # 1: recombine mesh: hex ;  0: no: only prisms
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('np(i)         : '+str(nb_plies)+'\n')     # 6, 12 or 24 # TODO: find a clever way of setting stacking sequence
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
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~MESH~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('lc(f)    : 1\n')        # mesh caracteristic size
	parameters.write('dx(i)    : 20\n')       # discretization in X direction
	parameters.write('ddy(i)   : 1\n')        # discretization of ply thickness
	parameters.write('dz(i)    : 60\n')       # discretization in Z direction
	parameters.write('dc(i)    : 7\n')        # discretization of corners
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~TRANSFORMATION~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('wrinkle(i)        : '+str(nb_wrinkles)+'\n') # Do we need to add wrinkle in the gridMod (c++) code, 2+ for multiple wrinkles
	parameters.write('Ramp(b)           : 1\n') #
	parameters.write('Abaqus_output(b)  : 0\n') #
	parameters.write('Dune_output(b)    : 1\n') #
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~WRINKLES~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('WID(s)             : Defect_'+str(cnt)+'\n')
	parameters.write('Wsize(f)    : '+str(p1)+'\n') # Amplitude max
	parameters.write('Wpos(f)     : 84.0,73.3,207.7\n') # center
	parameters.write('Wori(f)     : '+str(p4)+'\n') # Orientation in degree
	parameters.write('Wdamp(f)    : 8.0,'+str(p2)+','+str(p3)+'\n') # reduction of the amplitude through each direction
	# for i in range(nb_wrinkles):
	# 	if i<10:
	# 		parameters.write('Wsize'+str(i)+'(f)    : '+str(random.uniform(minWsize, maxWsize))+'\n') # Amplitude max
	# 		parameters.write('Wpos'+str(i)+'(f)     : '+str(random.uniform(minWposX, maxWposX))+','+str(random.uniform(minWposY, maxWposY))+','+str(random.uniform(minWposZ, maxWposZ))+'\n') # center
	# 		parameters.write('Wori'+str(i)+'(f)     : '+str(random.uniform(minWori, maxWori))+'\n') # Orientation in degree
	# 		parameters.write('Wdamp'+str(i)+'(f)    : '+str(random.uniform(minWdampX, maxWdampX))+','+str(random.uniform(minWdampY, maxWdampY))+','+str(random.uniform(minWdampZ, maxWdampZ))+'\n') # reduction of the amplitude through each direction
	# 	else:
	# 		parameters.write('Wsize'+str(i)+'(f)   : '+str(random.uniform(minWsize, maxWsize))+'\n') # Amplitude max
	# 		parameters.write('Wpos'+str(i)+'(f)    : '+str(random.uniform(minWposX, maxWposX))+','+str(random.uniform(minWposY, maxWposY))+','+str(random.uniform(minWposZ, maxWposZ))+'\n') # center
	# 		parameters.write('Wori'+str(i)+'(f)    : '+str(random.uniform(minWori, maxWori))+'\n') # Orientation in degree
	# 		parameters.write('Wdamp'+str(i)+'(f)   : '+str(random.uniform(minWdampX, maxWdampX))+','+str(random.uniform(minWdampY, maxWdampY))+','+str(random.uniform(minWdampZ, maxWdampZ))+'\n') # reduction of the amplitude through each direction
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~RAMP~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Rsize(f)          : 6.25\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~ABAQUS~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Path2result(s)    : Abaqus/results/\n')
	parameters.write('AbaqusOdbName(s)  : model\n')


	parameters.close()

def plot_S_dampY(x, title, colors):
    fig, ax = plt.subplots()
    for datax,datay,color in zip(np.array(x)[:, 0], np.array(x)[:, 1],colors):
        plt.scatter(datax,datay,color=color)
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
    for datax,datay,color in zip(np.array(x)[:, 0], np.array(x)[:, 2],colors):
        plt.scatter(datax,datay,color=color)
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
    for datax,datay,color in zip(np.array(x)[:, 0], np.array(x)[:, 3],colors):
        plt.scatter(datax,datay,color=color)
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
    for datax,datay,color in zip(np.array(x)[:, 1], np.array(x)[:, 2],colors):
        plt.scatter(datax,datay,color=color)
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

	n_samples = 100
	space = Space([(2., 12.), (1., 6.), (3., 10.), (-5., 5.)])

	lhs = Lhs(lhs_type="classic", criterion=None)
	x = lhs.generate(space.dimensions, n_samples)

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

	for cnt in range(n_samples):
		write_parameters(cnt,np.array(x)[cnt, 0], np.array(x)[cnt, 1], np.array(x)[cnt, 2], np.array(x)[cnt, 3])
