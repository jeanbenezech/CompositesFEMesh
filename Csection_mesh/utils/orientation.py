import numpy as np

from utils.parameters import *

class gridMod():
	def __init__(self,
				 u=None,
				 v=None,
				 label=None):
		self.u=u
		self.v=v
		self.label=label
	
	def init(self):
		self.u = np.zeros(3)
		self.v = np.zeros(3)
		self.label = 0
	
	def set_label(self, label):
		self.label=label
