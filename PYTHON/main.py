import numpy as np

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

	geo.dump(param)
