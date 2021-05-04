import numpy as npy
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('size', metavar='n', type=int, default=60,
                    help='size of the square mesh (default 60)')
parser.add_argument('fname', metavar='f', type=str)

args = parser.parse_args()
grid_size = (args.size, args.size)

m = npy.load('output-fdtd-2d-{}.npy'.format(args.fname))
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
