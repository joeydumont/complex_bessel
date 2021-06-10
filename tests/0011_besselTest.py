import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

file = h5py.File("hankelOverflow.h5", "r")

z = np.array(file["Coordinates"]["z"])
errorCode = np.array(file["Coordinates"]["errorCode"])
hankelH2 = np.array(file["BesselFunctionValues"]["hankelH2"], dtype=np.complex)
sumBessel = np.array(file["BesselFunctionValues"]["sumBessel"], dtype=np.complex)

diff = np.transpose(np.abs(hankelH2 - sumBessel))

print(errorCode)

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
for i in range(0,5):#range(np.size(sumBessel, 0)):
    ax.plot(z, diff[:, i])
    ax.set_ylim(0,1)
    ax2.plot(z, errorCode[i,:], label='{}'.format(i))
    plt.legend()

# ax.set_yscale('log')
ax.set_xscale("log")
plt.show()

for i in range(np.size(sumBessel, 0)):
  if max(diff[:, i]) > 0.5:
    print(i)
    break
