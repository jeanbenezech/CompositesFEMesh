import numpy as np
np.random.seed(123)

import sys
import random
import matplotlib.pyplot as plt
from skopt.space import Space
from skopt.sampler import Lhs
import matplotlib.cm as cm

# geometry in millimeters

def write_parameters(cnt=-1,p1=0, p2=0, p3=0, p4=0):
	nb_plies = 8
	nb_wrinkles = 1

	basename = 'w'
	i = cnt
	# for i in range(200):
	name = basename+str(i)

	parameters=open('wrinkles_parameters_'+str(i)+'.txt','w')
	parameters.write('Name: '+name+'\n') # Amplitude max
	parameters.write('Part_type: 0\n')
	parameters.write('Wrinkle_type: 0\n')
	parameters.write('Coord: '+str(47.0)+','+str(p1)+','+str(p2)+'\n') # center
	parameters.write('Orient: '+str(0.0)+'\n') # Orientation in degree
	parameters.write('Size: '+str(1.7)+'\n') # Orientation in degree
	parameters.write('Damp: '+str(0.18)+','+str(0.03)+','+str(0.09)+'\n') # reduction of the amplitude through each direction
	parameters.close()

def plot_S_dampY(x, title, colors):
    fig, ax = plt.subplots()
    for datax,datay,color in zip(np.array(x)[:, 0], np.array(x)[:, 1],colors):
        plt.scatter(datax,datay,color=color)
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 1], 'bo', label='samples')
    # plt.plot(np.array(x)[:, 0], np.array(x)[:, 1], 'bo', markersize=80, alpha=0.5)
    # ax.legend(loc="best", numpoints=1)
    ax.set_xlabel("Y")
    ax.set_xlim([-20., 160.])
    ax.set_ylabel("Z")
    ax.set_ylim([120., 360.])
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

	n_samples = 10
	space = Space([(160., 180.), (120., 360.)])
	lhs = Lhs(lhs_type="classic", criterion=None)
	x = lhs.generate(space.dimensions, n_samples)

	colors = cm.coolwarm(np.linspace(0, 1, n_samples))
	plot_S_dampY(x, 'Y vs Z',colors)
	# plot_S_dampZ(x, 'Size vs dampZ',colors)
	# plot_S_Angle(x, 'Size vs Angle',colors)
	# plot_DampYZ(x, 'dampY vs dampZ',colors)

	f=open('LhsWrinkles_2.txt', 'w')
	f.write('Y, Z\n')
	for i in range(n_samples):
		f.write(str(np.array(x)[i, 0])+', '+str(np.array(x)[i, 1])+'\n')
		# f.write(str(np.array(x)[i, 0])+', '+str(np.array(x)[i, 1])+', '+str(np.array(x)[i, 2])+', '+str(np.array(x)[i, 3])+'\n')
	f.close()

	for cnt in range(n_samples):
		write_parameters(cnt+100,np.array(x)[cnt, 0], np.array(x)[cnt, 1])


	# write_parameters()