# Python script to generate the contour plot
# seen in Abramowitz & Stegun's book on p. 359. 
# The values are imported from the file "contours.dat"
#
# The pylab module is required for this script to run
#
# Joey Dumont <joey.dumont@gmail.com>
# Denis Gagnon <gagnon88@gmail.com>
#

from pylab import *

Data = loadtxt('contours.dat')
x = Data[:,0]
y = Data[:,1]
M = sqrt(Data[:,2]**2 + Data[:,3]**2)

Z = Data[:,2]+complex(0,1)*Data[:,3]

phi=(180/pi)*abs(arctan2(Data[:,3],Data[:,2]))

Dimension = sqrt(M.size)

X=linspace(x.min(),x.max(), Dimension)
Y=linspace(y.min(),y.max(), Dimension)

Xg,Yg=meshgrid(X,Y)

M0=reshape(M,[Dimension, Dimension])
phi0=reshape(phi,[Dimension, Dimension])

contourM=linspace(0.2,3.2,16)
contourP=linspace(0,360,15)

figure(figsize=(7,5))

contour(Xg,Yg,M0,contourM)
CS = contour(Xg,Yg,phi0,contourP,colors='k',linestyles='dashdot')

Xcut=[-4.0,0]
Ycut=[0,0]

plot(Xcut,Ycut,lw=2.5,color='k')

xlabel('x')
ylabel('y')
title('Contour lines of the modulus and phase of $H_0^{(1)}(x+iy)$ \n (reproduced from Abramowitz \& Stegun, p.359)')

savefig('contours.png')