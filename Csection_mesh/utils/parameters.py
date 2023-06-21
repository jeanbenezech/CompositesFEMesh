import numpy as np

class sub():
	def __init__(self,
				 L=None, # lenght
				 thick=None, # local thickness
				 Ori=None, # gridMod (straight) max angle (curved)
				 DX=None,
				 shape=None):
		self.L = 0
		self.thick = 0
		self.Ori = 0
		self.DX = 0

class parameters():
	def __init__(self,
				 AbaqusName=None,
				 Shape=None,
				 resin=None,
				 resin_from_StackSeq=None,
				 CZ=None,
				 recombine=None,
				 nbp=None,
				 name=None,
				 dx=None,
				 dy=None,
				 ddy=None,
				 dflange=None,
				 dz=None,
				 dz_intervals=None,
				 dc=None,
				 nbp_by_surf=None,
				 nbl_by_surf=None,
				 X=None,
				 Height=None,
				 Y=None,
				 R=None,
				 Z=None,
				 e=None,
				 fl=None,
				 theta=None,
				 lc=None,
				 refOriPoint=None,
				 substr=None,
				 StackSeq=None,
				 ramp=None):
		self.nbp = nbp
		self.theta=0

	def init(self, filename):
		read_parameters(self, filename)

		self.refOriPoint= []
		self.substr = []

		substr = sub()
		substr.L = self.Height
		substr.T = 0
		substr.Ori = 0 * np.pi/180.
		substr.DX = self.dflange
		substr.shape='straight'

		substr2 = sub()
		substr2.L = self.Y+self.R
		substr2.T = 0
		substr2.Ori = 90 * np.pi/180.
		substr2.shape='curved'

		substr4 = sub()
		substr4.L = self.X
		substr4.T = 0
		substr4.Ori = 0 * np.pi/180.
		substr4.DX = self.dx
		substr4.shape='straight'

		if self.Shape==0: # Csection
			self.substr.append(substr)
			self.substr.append(substr2)
			self.substr.append(substr4)
			self.substr.append(substr2)
			self.substr.append(substr)
		elif self.Shape==1: # Flat laminate
			self.substr.append(substr4)

		if self.resin and self.CZ:
			self.dy = 4*self.nbp-2 # 1 ply = 1 comp layer + 2 resin elastic layer + 1 CZ (-3 at the end because no resin at the top)
		elif self.resin or self.CZ:
			self.dy = 2*self.nbp
		else:
			self.dy = (self.nbp+1)

	def write_input_orientation(self, filename):
		f=open(filename, 'w')

		f.write(str(len(self.refOriPoint))+'\n')
		for point in self.refOriPoint:
			f.write(str(point[0])+' '+str(point[1])+'\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read(filename):
	parameters = {}
	f=open(filename+'.txt', 'r')
	line=f.readline()
	while line!='':
		if line[0]!='~':
			parameters["{}".format(line.split('(')[0])] = "{}".format(line.strip().split()[2])
		line=f.readline()
	f.close()
	return parameters

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_parameters(param, filename):
	dico = read(filename)

	param.Shape = int(dico['Shape'])
	param.resin = int(dico['Auto_Resin_betw_plies'])
	param.recombine = int(dico['recombine'])
	param.CZ = int(dico['cohezive_elements'])
	# if param.resin==0:
	# 	param.CZ = 0

	param.nbp = int(dico['np'])
	param.X = float(dico['X'])
	param.Height = float(dico['Height'])
	param.Y = float(dico['Y'])
	param.R = float(dico['R'])
	param.Z = float(dico['Z'])
	param.e = float(dico['e'])

	param.lc = float(dico['lc'])

	param.name = dico['name']

	param.AbaqusName = dico['AbaqusOdbName']

	param.dx = int(dico['dx'])
	param.ddy = int(dico['ddy'])
	param.dflange = int(dico['dflange'])
	param.dc = int(dico['dc'])

	string_dz = dico['dz']
	if ',' not in string_dz:
		param.dz = np.asarray([int(dico['dz'])])
	else:
		param.dz = []
		print (string_dz.strip().split(','))
		for string_value in string_dz.strip().split(','):
			param.dz.append(int(string_value))

	param.ramp = int(dico['Ramp'])
	if param.ramp:
		string_StartEndinZdir = dico['StartEndinZdir']
		val = []
		for string_value in string_StartEndinZdir.strip().split(','):
			val.append(float(string_value))

		# print (val)

		param.dz_intervals = []
		interval = val[0]/float(param.Z) # start of the ramp in term of perone of the ful length
		param.dz_intervals.append(interval)
		interval += (2*val[1]+val[2])/float(param.Z) # length of the joggle region
		param.dz_intervals.append(interval)
		param.dz_intervals.append(1.)


	param.StackSeq = []
	for i in range(param.nbp):
		lplies=dico['p'+str(i)]
		tmp = np.asarray([float(lplies.strip().split(',')[0]),float(lplies.strip().split(',')[1]),float(lplies.strip().split(',')[2])])
		param.StackSeq.append(tmp)

	return param
