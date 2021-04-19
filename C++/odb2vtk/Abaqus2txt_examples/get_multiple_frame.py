import numpy as np
import sys
import os

# frames = np.asarray([])
frames = np.asarray(range(1, 8))

command = 'mkdir all_result'
os.system(command)

for f in frames:
    command = 'abaqus python GetCellFields.py model.odb '+ str(f)
    os.system(command)
    command = 'abaqus python Getdisplacement.py model.odb '+ str(f)
    os.system(command)
    command = 'mv model**.txt all_result/'
    os.system(command)
