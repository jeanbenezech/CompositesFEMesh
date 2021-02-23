import numpy as np

from utils.geometry import *
from utils.parameters import *
from utils.lines import *

def front_surface(param, geo):

	rep = reference()
	geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v = Clines(param, geo, rep)

	return geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v