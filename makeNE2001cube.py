from ne2001 import density_2001
import numpy as np
import sys
import pickle
import subprocess

if not len(sys.argv) == 5:
    print "Usage: step_x step_y step_z output"
    exit()

dx,dy,dz = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
x = np.linspace(-17.4,17.4,dx)
y = np.linspace(-17.4,17.4,dy)
z = np.linspace(0.,19.4,dz)

ne = np.zeros((dx,dy,2. * dz), dtype = np.float32)

for i in range(dx):
    for j in range(dy):
	for k in range(dz):
	    n = density_2001(x[i],y[j],z[k])
	    ne[i,j,k + dz - 1] =  n
	    if k:
		ne[i,j,dz - 1 - k] =  n

f = open(sys.argv[-1],'w')
pickle.dump(ne,f)
f.close()

subprocess.call(['gzip','-f',sys.argv[-1]])

#import matplotlib.pyplot as plt

#print ne.shape
#im = plt.imshow(np.log10(ne[:,:,dz]), aspect = 'auto')

#plt.colorbar(im)
#plt.show()
