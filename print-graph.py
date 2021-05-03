import numpy as npy
import matplotlib.pyplot as plt

grid_size = (60, 60)

m = npy.load('output-2d-s.npy')
m = m.reshape(-1, *grid_size)

x = abs(m).max()

plt.ion()
plt.show()

for t, grid in enumerate(m):
    plt.clf()
    plt.title("Timestep "+str(t))
    plt.imshow(grid, aspect='auto', cmap='seismic')
    plt.colorbar()
    plt.clim([-x, x])
    plt.draw()
    plt.pause(0.2)
