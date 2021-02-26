import numpy as np

class parameter():
	def __init__(self,
				 Lshape=None,
				 flat_limb=None,
				 resin=None,
				 CZ=None,
				 recombine=None,
				 is_im=None,
				 nbp=None,
				 nbp_c=None,
				 name=None,
				 dx=None,
				 dy=None,
				 dy_c=None,
				 ddy=None,
				 ddy_r=None,
				 dz=None,
				 ddz=None,
				 dc=None,
				 dl=None,
				 dside=None,
				 nbp_by_surf=None,
				 nbl_by_surf=None,
				 nb_total_plies=None,
				 reso=None,
				 X=None,
				 Y=None,
				 R=None,
				 Z=None,
				 tX=None,
				 tY=None,
				 tZ=None,
				 e=None,
				 fl=None,
				 limb=None,
				 side=None,
				 lc=None,
				 theta=None):
		self.nbp = nbp

	def init(self):
		read_parameters(self)

		if self.resin:
			if self.CZ: # 1 ply = 1 comp layer + 2 resin elastic layer + 1 CZ (-1 at the end because the filler is a resin layer)
				self.dy = 4*self.nbp
				self.dy_c = 4*self.nbp_c
			else: # We want a resin layer at the top (+1)
				self.dy = 2*self.nbp+1
				self.dy_c = 2*self.nbp_c+1
		else:
			self.dy = (self.nbp+1)
			self.dy_c = (self.nbp_c+1)

		# To define the size of Volume vector
		if(self.flat_limb==0):
			self.nb_total_plies = (6*(self.dy_c-1)+3*(self.dy-1)) * self.dz + 4 # 3* three laminate (bottom left and right) and +4 for the filler and each limits of the filler
		else:
			self.nb_total_plies = (6*(self.dy_c-1)+5*(self.dy-1)) * self.dz + 4 # 3* laminate (bottom, up left and up right) and 2 flat sides // and +4 for the filler and each limits of the filler

	def write_input_orientation(self, filename):
		f=open(filename, 'w')
		val = 1
		P0 = np.asarray([0.0, 0.0])
		P1 = np.asarray([0.0, 0.0])

		P0 = np.asarray([self.limb, self.Y])
		if self.flat_limb:
			P0[0]+=self.side
		if self.CZ:
			P0[1]-=self.e

		P1 = np.asarray([P0[0] + 2*self.X + self.fl, P0[1] + self.X + self.fl])
		if self.CZ:
			P1[1]+=self.e

		typ = 1
		if (self.Lshape==1):
			typ = 4
		
		f.write(str(val)+' '+str(typ)+' '+str(self.theta)+'\n')
		f.write(str(P0[0])+' '+str(P0[1])+'\n')
		f.write(str(P1[0])+' '+str(P1[1])+'\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read():
	parameters = {}
	f=open('parameters.txt', 'r')
	line=f.readline()
	while line!='':
		if line[0]!='~':
			parameters["{}".format(line.split('(')[0])] = "{}".format(line.strip().split()[2])
		line=f.readline()
	f.close()
	return parameters

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_parameters(param):
	dico = read()

	param.Lshape = int(dico['L-shape'])
	param.flat_limb = int(dico['Flat_Limb'])
	param.resin = int(dico['Resin_bet_plies'])
	param.recombine = int(dico['recombine'])
	param.CZ = int(dico['cohezive_elements'])
	if param.resin==0:
		param.CZ = 0

	param.nbp = int(dico['np'])
	param.nbp_c = int(dico['np_c'])
	param.X = float(dico['X'])
	param.Y = float(dico['Y'])
	param.R = float(dico['R'])
	param.Z = float(dico['Z'])
	param.e = float(dico['e'])
	param.fl = float(dico['fl'])
	param.limb = float(dico['limb'])

	param.side = float(dico['side'])

	param.lc = float(dico['lc'])

	param.name = dico['name']

	param.dx = int(dico['dx'])
	param.ddy = int(dico['ddy'])
	param.ddy_r = int(dico['ddy_r'])
	param.dz = int(dico['dz'])
	param.ddz = int(dico['ddz'])
	param.dc = int(dico['dc'])
	param.dl = int(dico['dl'])
	param.dside = int(dico['dside'])

	param.theta = float(dico['hat_angle'])
	if (param.Lshape==0):
		param.theta=90.0

	return param
