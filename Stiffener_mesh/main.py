import numpy as np
import os

from utils.surfaces import *
from utils.parameter import *
from utils.geometry import *

if __name__ == '__main__':

	param = parameter()
	param.init()
	param.write_input_orientation('input.txt')

	geo = geometry()
	geo.init(param)

	geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v, N = front_surface(param, geo)

	if param.Lshape==1:
		filler_L_shape(param, geo, N)
	else:
		filler(param, geo, N)

	# Ajust coordinate system
	# geo.ajust(param)

	geo.dump(param)

	command = 'gmsh {} -3 -format msh2 -o {}'.format(param.name+'.geo', param.name+'.msh')
	os.system(command)

