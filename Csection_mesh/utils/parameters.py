import numpy as np
import configparser
import copy

class sub():
	def __init__(self,
				 L=None, # lenght
				 thick=None, # local thickness
				 Ori=None, # gridMod (straight) max angle (curved)
				 DX=None,
				 FlangeRotAngle=None,
				#  FlangeBC = None,
				 shape=None,
				 type=None):
		self.L = 0
		self.Linit = 0
		self.A = 0
		self.thick = 0
		self.Ori = 0
		self.DX = 0
		# self.FlangeBC = 0
		self.FlangeRotAngle = 0
		self.type = 'classic'

	def local_Ix0(self, param, d):
		local_xmax = param.Height
		local_radius = param.R + d
		theta = -np.pi/2. + (self.FlangeRotAngle * np.pi/180.)
		Cx = local_xmax
		# self.A=0
		Ix = Cx - (( Cx + local_radius * np.cos(theta)) / np.sin(-theta))
		return Ix

	def local_Ix(self, param, d):
		local_xmax = param.Height
		local_radius = param.R + d
		theta = -np.pi/2. + (self.FlangeRotAngle * np.pi/180.)
		Cx = local_xmax
		Ix = Cx - (( Cx + local_radius * np.cos(theta) - self.A) / np.sin(-theta))
		return Ix


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
				 RotateFlanges=None,
				 AngleRotateFlangeRi=None,
				 AngleRotateFlangeLe=None,
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

			substr_flange_BC = []

			substr_flange_BC1 = sub()
			substr_flange_BC1.type = 'flangeRotLe'
			substr_flange_BC1.FlangeRotAngle = self.AngleRotateFlangeLe
			substr_flange_BC1.Linit = self.dflange_intervals[0]
			substr_flange_BC1.A = self.dflange_intervals[0]
			substr_flange_BC1.T = 0
			substr_flange_BC1.Ori = 0 * np.pi/180.
			substr_flange_BC1.DX = self.dflangeBot
			substr_flange_BC1.shape='straight'
			substr_flange_BC.append(substr_flange_BC1)

			for iter_BC in range(1,len(self.dflange_intervals)):
				substr_flange_BCi = sub()
				substr_flange_BCi.type = 'flangeRotLe'
				substr_flange_BCi.FlangeRotAngle = self.AngleRotateFlangeLe
				substr_flange_BCi.Linit =  self.dflange_intervals[iter_BC]-self.dflange_intervals[iter_BC-1]
				substr_flange_BCi.A = self.dflange_intervals[iter_BC]
				substr_flange_BCi.T = 0
				substr_flange_BCi.Ori = 0 * np.pi/180.
				substr_flange_BCi.DX = 2
				substr_flange_BCi.shape='straight'
				substr_flange_BC.append(substr_flange_BCi)

			substr_flange_end = sub()
			substr_flange_end.L = self.Height - self.dflange_intervals[-1] #+ y
			substr_flange_end.T = 0
			substr_flange_end.Ori = 0 * np.pi/180.
			substr_flange_end.DX = self.dflangeTop
			substr_flange_end.shape='straight'

			substr_flange_down_BC = []

			substr_flange_beg = sub()
			substr_flange_beg.type = 'flangeRotRi'
			substr_flange_beg.FlangeRotAngle = self.AngleRotateFlangeRi
			substr_flange_beg.Linit = self.Height - self.dflange_intervals[-1]
			substr_flange_beg.A = self.dflange_intervals[-1]
			substr_flange_beg.T = 0
			substr_flange_beg.Ori = 0 * np.pi/180.
			substr_flange_beg.DX = self.dflangeTop
			substr_flange_beg.shape='straight'
			substr_flange_down_BC.append(substr_flange_beg)

			for iter_BC in range(1,len(self.dflange_intervals)):
				substr_flange_BCi = sub()
				substr_flange_BCi.type = 'flangeRotRi'
				substr_flange_BCi.FlangeRotAngle = self.AngleRotateFlangeRi
				substr_flange_BCi.Linit = - self.dflange_intervals[-iter_BC-1] + self.dflange_intervals[-iter_BC]
				substr_flange_BCi.A = self.dflange_intervals[-iter_BC-1]
				substr_flange_BCi.T = 0
				substr_flange_BCi.Ori = 0 * np.pi/180.
				substr_flange_BCi.DX = 2
				substr_flange_BCi.shape='straight'
				substr_flange_down_BC.append(substr_flange_BCi)

			substr_flange_down_end = sub()
			substr_flange_down_end.type = 'flangeRotRi'
			substr_flange_down_end.FlangeRotAngle = self.AngleRotateFlangeRi
			substr_flange_down_end.Linit = self.dflange_intervals[0]
			substr_flange_down_end.A = 0.
			substr_flange_down_end.T = 0
			substr_flange_down_end.Ori = 0 * np.pi/180.
			substr_flange_down_end.DX = self.dflangeBot
			substr_flange_down_end.shape='straight'
			# # substr_flange_down_end.L = self.dflange_intervals[0]
			# substr_flange_down_end.T = 0
			# substr_flange_down_end.Ori = 0 * np.pi/180.
			# substr_flange_down_end.shape='straight'

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
			for iter_BC in range(len(substr_flange_BC)):
				self.substr.append(substr_flange_BC[iter_BC])
			self.substr.append(substr_flange_end)
			self.substr.append(substr2)
			self.substr.append(substr4)
			self.substr.append(substr2)
			for iter_BC in range(len(substr_flange_BC)):
				self.substr.append(substr_flange_down_BC[iter_BC])
			self.substr.append(substr_flange_down_end)
			# self.substr.append(substr_flange_1)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_3)
			# self.substr.append(substr2)
			# self.substr.append(substr4)
			# self.substr.append(substr2)
			# self.substr.append(substr_flange_3)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_2)
			# self.substr.append(substr_flange_1)
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

	param.RotateFlanges = config['GEO-TRANSFORMATION'].getint('RotateFlanges(b)')
	param.AngleRotateFlangeRi = config['GEO-TRANSFORMATION'].getfloat('AngleRotateFlangeRi(f)')
	param.AngleRotateFlangeLe = config['GEO-TRANSFORMATION'].getfloat('AngleRotateFlangeLe(f)')
	if(param.RotateFlanges==0):
		param.AngleRotateFlangeRi= 0.
		param.AngleRotateFlangeLe= 0.

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

		# print (val)

	
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
