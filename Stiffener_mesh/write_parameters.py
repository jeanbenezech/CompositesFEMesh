import sys

# geometry in millimeters

def write_parameters():

	parameters =open('parameters.txt','w')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('name(s)                : stiffener\n') # mesh name
	parameters.write('L-shape(b)             : 0\n')   # 1: stiffener L-shape ; 0 stiffener T-shape
	parameters.write('Flat_Limb(b)           : 1\n')   # 1: yes ; 0: no
	parameters.write('Resin_bet_plies(b)     : 1\n')   # 1: yes ; 0: no
	parameters.write('cohezive_elements(b)   : 0\n')   # 1: yes ; 0: no
	parameters.write('recombine(b)           : 1\n')   # 1: recombine mesh: hex, prisms ;  0: no: only prisms
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('np_c(i)       : 5\n')       # Number of plies by curved laminate
	parameters.write('np(i)         : 4\n')       # Number of plies for the flat plate
	parameters.write('X(f)          : 6.0\n')    #
	parameters.write('Y(f)          : 4.0\n')    #
	parameters.write('R(f)          : 2.0\n')    #
	parameters.write('Z(f)          : 3.0\n')    #
	parameters.write('e(f)          : 0.01\n')  # Resin layer thickness
	parameters.write('fl(f)         : 0.05\n')  # Filler limbs thickness
	parameters.write('limb(f)       : 2.0\n')    #
	parameters.write('side(f)       : 3.0\n')    #
	parameters.write('hat_angle(f)  : 90.0\n')   #
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~MESH~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('lc(f)    : 1.0\n') # mesh caracteristic size
	parameters.write('dx(i)    : 15\n')  #
	parameters.write('ddy(i)   : 2\n')   #
	parameters.write('ddy_r(i) : 2\n')   #
	parameters.write('dz(i)    : 1\n')   #
	parameters.write('ddz(i)   : 4\n')   #
	parameters.write('dc(i)    : 12\n')  #
	parameters.write('dl(i)    : 5\n')   #
	parameters.write('dside(i) : 5\n')   #
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~TRANSFORMATION~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('wrinkle(b)        : 1\n') # Do we need to add wrinkle in the gridMod (c++) code
	parameters.write('Ramp(b)           : 0\n') #
	parameters.write('Abaqus_output(b)  : 1\n') #
	parameters.write('Dune_output(b)    : 0\n') #
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~WRINKLES~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Wsize(f)          : 4.0\n') # Amplitude max
	parameters.write('Wpos(f)           : 86.0,75.0,250.0\n') # center
	parameters.write('Wori(f)           : 30.0\n') # Orientation in degree
	parameters.write('Wdamp(f)          : 2.0,0.7,10.0\n') # reduction of the amplitude through each direction
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
