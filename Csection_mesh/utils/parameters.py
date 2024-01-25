import numpy as np
import configparser

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
				 dflangeBot=None,
				 dflangeTop=None,
				 dz=None,
				 dz_intervals=None,
				 dflange_intervals=None,
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
				 ramp=None,
				 RigidBoundary=None):
		self.nbp = nbp
		self.theta=0

	def init(self, filename):
		read_ini_parameters(self, filename)

		self.refOriPoint= []
		self.substr = []

		if self.Shape==0: # Csection
			substr = sub()
			substr.L = self.Height
			substr.T = 0
			substr.Ori = 0 * np.pi/180.
			substr.DX = self.dflange
			substr.shape='straight'

		##### Sub - substr of the flange for the Specific BC
		if self.Shape==2: # Csection with specified position of node for BC
			substr_flange_1 = sub()
			substr_flange_1.L = self.dflange_intervals[0]
			substr_flange_1.T = 0
			substr_flange_1.Ori = 0 * np.pi/180.
			substr_flange_1.DX = self.dflangeBot
			substr_flange_1.shape='straight'

			substr_flange_2 = sub()
			substr_flange_2.L = 2
			substr_flange_2.T = 0
			substr_flange_2.Ori = 0 * np.pi/180.
			substr_flange_2.DX = 2
			substr_flange_2.shape='straight'

			substr_flange_3 = sub()
			substr_flange_3.L = self.Height - self.dflange_intervals[-1]
			substr_flange_3.T = 0
			substr_flange_3.Ori = 0 * np.pi/180.
			substr_flange_3.DX = self.dflangeTop
			substr_flange_3.shape='straight'

		if self.Shape==2 or self.Shape==0: # Csection with specified position of node for BC
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
		elif self.Shape==2: # Csection with specified position of node for BC
			# self.substr.append(substr)
			self.substr.append(substr_flange_1)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_3)
			self.substr.append(substr2)
			self.substr.append(substr4)
			self.substr.append(substr2)
			self.substr.append(substr_flange_3)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_2)
			self.substr.append(substr_flange_1)
			# self.substr.append(substr)

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
def read_ini_parameters(param, filename): #### deprecated
	# dico = read(filename)
	config = configparser.ConfigParser()
	config.read(filename+'.ini')

	param.Shape = config['GENERAL'].getint('Shape(i)')
	param.resin = config['GENERAL'].getint('Auto_Resin_betw_plies(b)')
	param.recombine = config['GENERAL'].getint('recombine(b)')
	param.CZ = config['GENERAL'].getint('cohezive_elements(b)')
	# if param.resin==0:
	# 	param.CZ = 0

	param.nbp = config['GEOMETRY'].getint('np(i)')
	param.X = config['GEOMETRY'].getfloat('X(f)')
	param.Height = config['GEOMETRY'].getfloat('Height(f)')
	param.Y = config['GEOMETRY'].getfloat('Y(f)')
	param.R = config['GEOMETRY'].getfloat('R(f)')
	param.Z = config['GEOMETRY'].getfloat('Z(f)')
	param.e = config['GEOMETRY'].getfloat('e(f)')

	param.lc = config['MESH'].getfloat('lc(f)')

	param.name = config['GENERAL']['name(s)']

	param.AbaqusName = config['GENERAL']['AbaqusOdbName(s)']

	param.dx = config['MESH'].getint('dx(i)')
	param.ddy = config['MESH'].getint('ddy(i)')
	param.dflange = config['MESH'].getint('dflange(i)')
	param.dflangeBot = config['MESH'].getint('dflange_bot(i)')
	param.dflangeTop = config['MESH'].getint('dflange_top(i)')
	param.dc = config['MESH'].getint('dc(i)')
	string_dz = config['MESH']['dz(i)']
	if ',' not in string_dz:
		param.dz = np.asarray([int(string_dz)])
	else:
		param.dz = []
		print (string_dz.strip().split(','))
		for string_value in string_dz.strip().split(','):
			param.dz.append(int(string_value))

	param.ramp = config['GEO-TRANSFORMATION'].getint('Ramp(b)')
	if param.ramp:
		string_StartEndinZdir = config['GEO-TRANSFORMATION']['StartEndinZdir(f)']
		val = []
		for string_value in string_StartEndinZdir.strip().split(','):
			val.append(float(string_value))
		# print (val)

	param.RigidBoundary = config['GEO-TRANSFORMATION'].get('RigidBoundary(b)')
	string_dZInterval = config['GEO-TRANSFORMATION'].get('dZInterval(f)')
	if len(param.dz)>1 and param.RigidBoundary:
		param.dz_intervals = []
		for string_value in string_dZInterval.strip().split(','):
			val = float(string_value)
			interval = val/float(param.Z) # start of the ramp in term of percentage of the full length
			param.dz_intervals.append(interval)
		param.dz_intervals.append(1.)

	string_dFlangeInterval = config['GEO-TRANSFORMATION'].get('dFlangeInterval(b)')
	if param.Shape == 2:
		param.dflange_intervals = []
		for string_value in string_dFlangeInterval.strip().split(','):
			val = float(string_value)
			param.dflange_intervals.append(val)

	param.StackSeq = []
	for i in range(param.nbp):
		lplies=config['STACKSEQ'].get('p'+str(i)+'(f,f,b)')
		tmp = np.asarray([float(lplies.strip().split(',')[0]),float(lplies.strip().split(',')[1]),float(lplies.strip().split(',')[2])])
		param.StackSeq.append(tmp)

	return param

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read(filename): #### deprecated
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
def read_parameters(param, filename): #### deprecated
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
	param.dflangeBot = int(dico['dflange_bot'])
	param.dflangeTop = int(dico['dflange_top'])
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

		print (val)

	
	param.RigidBoundary = dico['RigidBoundary']
	string_dZInterval = dico['dZInterval']
	if len(param.dz)>1 and param.RigidBoundary:
		param.dz_intervals = []
		for string_value in string_dZInterval.strip().split(','):
			val = float(string_value)
			interval = val/float(param.Z) # start of the ramp in term of percentage of the full length
			param.dz_intervals.append(interval)
		param.dz_intervals.append(1.)

	string_dFlangeInterval = dico['dFlangeInterval']
	if param.Shape == 2:
		param.dflange_intervals = []
		for string_value in string_dFlangeInterval.strip().split(','):
			val = float(string_value)
			param.dflange_intervals.append(val)

	param.StackSeq = []
	for i in range(param.nbp):
		lplies=dico['p'+str(i)]
		tmp = np.asarray([float(lplies.strip().split(',')[0]),float(lplies.strip().split(',')[1]),float(lplies.strip().split(',')[2])])
		param.StackSeq.append(tmp)

	return param
