import os
import numpy as npy
import matplotlib.pyplot as plt

os.chdir('/Users/eric/hpc/projects/hpc-final/')
m = npy.load('output-2d-s.npy')
#plt.plot(m.T)
plt.imshow(m)
plt.colorbar()
plt.show()