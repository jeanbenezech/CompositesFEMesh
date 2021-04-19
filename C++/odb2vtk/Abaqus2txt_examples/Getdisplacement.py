import odbAccess  
from abaqusConstants import *
import numpy as np
import sys
import os

def export_field(filename, nframe=-1):
	
	odb = odbAccess.openOdb(filename)
	
	nvertex = len(odb.rootAssembly.instances['M'].nodes)
	
	U = odb.steps['Step-1'].frames[nframe].fieldOutputs['U'].values
	
	nameout = filename[:-4]+'_U.txt'
	if(nframe>=0):
		nameout = filename[:-4]+'_frame_'+str(nframe)+'_U.txt'
	
	g=open(nameout, 'w')
	g.write('nodeLabel, U1, U2, U3')
	g.write('\n')
	for u in U:
		g.write(str(u.nodeLabel))
		[g.write(' '+str(u.data[i])) for i in range(3)]
		g.write('\n')
	
if __name__ == '__main__':
	
	filename = sys.argv[1]
	
	nframe=-1
	if len(sys.argv)>2:
		nframe = int(sys.argv[2])
	print ('Get odb fields at frame number '+str(nframe))
	
	export_field(filename, nframe)
