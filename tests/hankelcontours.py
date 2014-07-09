# Python script to generate the contour plot
# seen in Abramowitz & Stegun's book on p. 359. 
# The values are imported from the file "contours.dat"
#
# The pylab module is required for this script to run
#
# Joey Dumont <joey.dumont@gmail.com>
# Denis Gagnon <gagnon88@gmail.com>
#

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 10
plt.rcParams['legend.numpoints'] = 3

Data = np.loadtxt('contours.dat')
x = Data[:,0]
y = Data[:,1]
M = np.sqrt(Data[:,2]**2 + Data[:,3]**2)

Z = Data[:,2]+np.complex(0,1)*Data[:,3]

phi=(180/np.pi)*np.abs(np.arctan2(Data[:,3],Data[:,2]))

Dimension = np.sqrt(M.size)

X=np.linspace(x.min(),x.max(), Dimension)
Y=np.linspace(y.min(),y.max(), Dimension)

Xg,Yg=np.meshgrid(X,Y)

M0=np.reshape(M,[Dimension, Dimension])
phi0=np.reshape(phi,[Dimension, Dimension])

contourM=np.linspace(0.2,3.2,16)
contourP=np.linspace(0,360,15)

plt.figure(figsize=(7,5))

plt.contour(Xg,Yg,M0,contourM)
CS = plt.contour(Xg,Yg,phi0,contourP,colors='k',linestyles='dashdot')

Xcut=[-4.0,0]
Ycut=[0,0]

plt.plot(Xcut,Ycut,lw=2.5,color='k')

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.title('Contour lines of the modulus and phase of $H_0^{(1)}(x+iy)$ \n (reproduced from Abramowitz \& Stegun, p.359)')

plt.savefig('contours.png')
