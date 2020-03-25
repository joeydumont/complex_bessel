import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import h5py

file = h5py.File("hankelOverflow.h5", 'r')

z = np.array(file['Coordinates']['z'])
hankelH2 = np.array(file['BesselFunctionValues']['hankelH2'], dtype=np.complex)
sumBessel = np.array(file['BesselFunctionValues']['sumBessel'], dtype=np.complex)

diff = np.transpose(np.abs(hankelH2-sumBessel))

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(np.size(sumBessel,0)):
  ax.plot(z,diff[:,i])

#ax.set_yscale('log')
ax.set_xscale('log')
plt.show()
