# -------------------- Information -------------------- #
# Author:       Joey Dumont <joey.dumont@gmail.com>     #
# Date created: Jun. 17th, 2014							#
# Date mod.:    Jun. 17th, 2014							#
# Description:  We plot tthe values of Ai(z) for large	#
#				arguments. 								#
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #
import matplotlib.pyplot as plt
import numpy as np
import morgenstemning as mrg
from matplotlib.colors import LogNorm
ms,msi = mrg.morgenstemning()

from scipy.special import airy

# Setting the rc parameters.
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}\usepackage[charter]{mathdesign}"]
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 10
plt.rcParams['legend.numpoints'] = 3

ticks = np.arange(0, 1200, 200)
labels = np.arange(-500, 700, 200)

nbPoints = 250
x = np.linspace(-500,500,2*nbPoints+1)
y = np.linspace(-500,500,2*nbPoints+1)
X, Y = np.meshgrid(x,y)

z = np.zeros((2*nbPoints+1,2*nbPoints+1),dtype=np.complex)

for i in range(2*nbPoints+1):
	for j in range(500):
		airyTest = airy(x[i]+1j*y[j])[0]
		if (airyTest == airyTest):
			z[i,j] = airy(x[i]+1j*y[j])[0]

figAiScipyRe = plt.figure(figsize=(4,4))
axAiScipyRe = figAiScipyRe.add_subplot(111)
plt.pcolormesh(X, Y, np.abs(np.real(z)), cmap=msi, norm=LogNorm(), rasterized=True)
plt.title(r"|Re$(\text{Ai}(x+iy))|$ -- SciPy")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
#ax1.xaxis.set_ticks(np.arange(-500,501,100))
#ax1.yaxis.set_ticks(np.arange(-500,501,100))
plt.savefig("realAiry-scipy.png", bbox_inches='tight')

figAiScipyIm = plt.figure(figsize=(4,4))
axAiScipyIm = figAiScipyIm.add_subplot(111)
plt.pcolormesh(abs(np.imag(z)), cmap=msi, norm=LogNorm(), rasterized=True)
plt.xticks(ticks,labels)
plt.yticks(ticks,labels)
plt.title(r"|Im$(\text{Ai}(x+iy))|$ -- SciPy")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
axAiScipyIm.set_xlim((0,1000))
axAiScipyIm.set_ylim((0,1000))
#ax1.xaxis.set_ticks(np.arange(-500,501,100))
#ax1.yaxis.set_ticks(np.arange(-500,501,100))
plt.savefig("imagAiry-scipy.png", bbox_inches='tight')



# ---------------- Data Importation ------------------- #
realAiry = np.loadtxt("realAiry.dat")
imagAiry = np.loadtxt("imagAiry.dat")
realAirynoOPT = np.loadtxt("realAiry-noOPT.dat")
imagAirynoOPT = np.loadtxt("imagAiry-noOPT.dat")

fig1 = plt.figure(figsize=(4,4))
ax1 = fig1.add_subplot(111)
plt.pcolormesh(abs(realAiry), cmap=msi, norm=LogNorm(), rasterized=True)
plt.xticks(ticks,labels)
plt.yticks(ticks,labels)
plt.title(r"|Re$(\text{Ai}(x+iy))|$ -- with -O3")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
ax1.set_xlim((0,1000))
ax1.set_ylim((0,1000))
#ax1.xaxis.set_ticks(np.arange(-500,501,100))
#ax1.yaxis.set_ticks(np.arange(-500,501,100))
plt.savefig("realAiry-opt.png", bbox_inches='tight')


fig2 = plt.figure(figsize=(4,4))
ax2 = fig2.add_subplot(111)
plt.pcolor(abs(realAirynoOPT), cmap=msi, norm=LogNorm(), rasterized=True)
plt.xticks(ticks,labels)
plt.yticks(ticks,labels)
plt.title(r"$|\text{Re}(\text{Ai}(x+iy))|$ -- with -O0")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
ax2.set_xlim((0,1000))
ax2.set_ylim((0,1000))
#ax1.xaxis.set_ticks(np.arange(-500,501,100))
#ax1.yaxis.set_ticks(np.arange(-500,501,100))
plt.savefig("realAiry-noOpt.png", bbox_inches='tight')

fig3 = plt.figure(figsize=(4,4))
ax3 = fig3.add_subplot(111)
plt.pcolor(abs(imagAiry), cmap=msi, norm=LogNorm(), rasterized=True)
plt.xticks(ticks,labels)
plt.yticks(ticks,labels)
plt.title(r"|Im$(\text{Ai}(x+iy))|$ -- with -O3")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
ax3.set_xlim((0,1000))
ax3.set_ylim((0,1000))
#ax1.xaxis.set_ticks(np.arange(-500,501,100))
#ax1.yaxis.set_ticks(np.arange(-500,501,100))
plt.savefig("imagAiry-opt.png", bbox_inches='tight')

fig4 = plt.figure(figsize=(4,4))
ax4 = fig4.add_subplot(111)
plt.pcolor(abs(imagAirynoOPT), cmap=msi, norm=LogNorm(), rasterized=True)
plt.xticks(ticks,labels)
plt.yticks(ticks,labels)
plt.title(r"|Im$(\text{Ai}(x+iy))|$ -- with -O0")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
ax4.set_xlim((0,1000))
ax4.set_ylim((0,1000))
#ax1.xaxis.set_ticks(np.arange(-500,501,100))
#ax1.yaxis.set_ticks(np.arange(-500,501,100))
plt.savefig("imagAiry-noOpt.png", bbox_inches='tight')
