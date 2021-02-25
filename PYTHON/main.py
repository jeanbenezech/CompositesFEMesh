import numpy as np
import os

from utils.surfaces import *
from utils.parameters import *
from utils.geometry import *

if __name__ == '__main__':

	param = parameters()
	param.init()

	geo = geometry()
	geo.init(param)

	geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v = front_surface(param, geo)

	param.write_input_orientation('input.txt')

	# Ajust coordinate system
	# geo.ajust(param)

	command = 'gmsh {} -3 -format msh2 -o {}'.format(param.name+'.geo', param.name+'.msh')
	os.system(command)

	geo.dump(param)
