import sys

# geometry in millimeters

def write_parameters():

	parameters=open('parameters.txt','w')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~GENERAL~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('Resin_betw_plies(b)    : 0\n')   # 1: yes ; 0: no
	parameters.write('cohezive_elements(b)   : 0\n')   # 1: yes ; 0: no
	parameters.write('recombine(b)           : 1\n')   # 1: recombine mesh: hex ;  0: no: only prisms
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~GEOMETRY~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('np(i)         : 24\n')    #
	parameters.write('X(f)          : 150\n')   #
	parameters.write('Y(f)          : 6.48\n')  #
	parameters.write('R(f)          : 10\n')    #
	parameters.write('Height(f)     : 75\n')    #
	parameters.write('Z(f)          : 500\n')   #
	parameters.write('e(f)          : 0.01\n')  # Resin layer thickness
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~MESH~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('name(s)  : Csection\n') # mesh name
	parameters.write('lc(f)    : 1\n')        # mesh caracteristic size
	parameters.write('dx(i)    : 12\n')       #
	parameters.write('ddy(i)   : 2\n')        #
	parameters.write('dz(i)    : 30\n')      #
	parameters.write('dc(i)    : 7\n')       #
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~TRANSFORMATION~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	parameters.write('wrinkle(b)        : 0\n') # Do we need to add wrinkle in the gridMod (c++) code
	parameters.write('Abaqus_output(b)  : 1\n') # Do we need to add wrinkle in the gridMod (c++) code
	parameters.write('Dune_output(b)    : 1\n') # Do we need to add wrinkle in the gridMod (c++) code

	parameters.close()

if __name__ == '__main__':

	write_parameters()
