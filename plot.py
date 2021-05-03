import matplotlib.pyplot as plt
import numpy as np


data = np.loadtxt("2d_data.csv", delimiter=",")
data = data.reshape(-1, 60, 60)
bound = max(abs(data.min()), abs(data.max()))

plt.ion()
plt.show()

for i, frame in enumerate(data):
    plt.clf()
    plt.imshow(frame, cmap='seismic')
    plt.colorbar()
    plt.clim([-bound, bound])
    plt.pause(0.001)
