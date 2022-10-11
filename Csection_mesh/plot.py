import numpy as np

import matplotlib.pyplot as plt

filename = "test.txt"

f = open(filename, "r")
lines = f.readlines()
X=[]
Y=[]
for l in lines:
    X.append(float(l.strip().split(' : ')[0]))
    Y.append(float(l.strip().split(' : ')[1]))


fig, ax = plt.subplots(figsize=(12, 9))
ax.plot(X, Y, 'o-',linewidth=3, markersize=10)
plt.savefig("position.png", dpi=120)