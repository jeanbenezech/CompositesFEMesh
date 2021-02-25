import odbAccess
from abaqusConstants import *
import numpy as np
import sys
import os

def export_field(filename, nframe=-1):

	odb = odbAccess.openOdb(filename)
	myCSYS=odb.rootAssembly.DatumCsysByThreePoints(name='myCSYS',
	coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0),
	point2=(0.0, 1.0, 0.0))

	#Export global values : S
	nelem = len(odb.rootAssembly.instances['M'].elements)

	field_global = np.zeros((nelem, 14))

	S = odb.steps['Step-1'].frames[nframe].fieldOutputs['S'].getTransformedField(myCSYS).getSubset(position=CENTROID).values
	E = odb.steps['Step-1'].frames[nframe].fieldOutputs['E'].getTransformedField(myCSYS).getSubset(position=CENTROID).values
	SDEG = odb.steps['Step-1'].frames[nframe].fieldOutputs['SDEG'].values
	QUAD = odb.steps['Step-1'].frames[nframe].fieldOutputs['QUADSCRT'].values
	
	for i in range(0, 6):
		for s in S:
			field_global[s.elementLabel-1, i] = s.data[i]
		for e in E:
			field_global[e.elementLabel-1, i+6] = e.data[i]
	for d in SDEG:
		field_global[d.elementLabel-1, 12] = d.data
		print(d.elementLabel)
	for q in QUAD:
		field_global[q.elementLabel-1, 13] = q.data

	nameout = filename[:-4]+'.txt'
	if(nframe>=0):
		nameout = filename[:-4]+'_frame_'+str(nframe)+'.txt'

	g=open(nameout, 'w')
	g.write('Nelem, S11, S22, S33, S12, S13, S23, E11, E22, E33, E12, E13, E23, SDEG, QUADSCRT')
	g.write('\n')
	for i,j in enumerate(field_global[:,0]):
		g.write(str(i))
		[g.write(' '+str(k)) for k in field_global[i,:]]
		g.write('\n')

if __name__ == '__main__':

	filename = sys.argv[1]

	nframe=-1
	if len(sys.argv)>2:
		nframe = int(sys.argv[2])
	print 'Get odb fields at frame number '+str(nframe)

	export_field(filename, nframe)
