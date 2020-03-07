# -------------------- Information -------------------- #
# Author:       Joey Dumont <joey.dumont@gmail.com>     #
# Created:      Nov. 26th, 2014                         #
# Modified:     Nov. 27th, 2014                         #
# Description:  We compute the values of the Bessel     #
#               and Hankel functions in a given range   #
#               and store them in a structured HDF5     #
#               file.                                   #
# Dependencies: - python>=3.2                           #
#               - OpenMPI                               #
#               - HDF5                                  #
#               - mpi4py                                #
#               - h5py (OpenMPI)                        #
#               - Scipy, NumPy                          #
# ----------------------------------------------------- #

# --------------- Modules Importation ----------------- #

# To exit early.
from sys import exit

# Math functions
import math
import numpy as np

# -- HDF5 with OpenMPI support
from mpi4py import MPI
import h5py

# -- Bessel functions.
from scipy.special import jvp, yvp, kvp, ivp
from scipy.special import h1vp, h2vp

# -- Argument Parser
import argparse

# --------------- Argument Processing ----------------- #

# -- We create the necessary arguments and parse them.
parser = argparse.ArgumentParser(prog='testBesselValues',
                                 description="Computes values of the Bessel functions for different orders and arguments.")

parser.add_argument("--filename", help="Name of the file to create.", required=True)
parser.add_argument("--pmax", type=int, default=0, help="Maximum order of differentiation of Bessel functions considered.")
parser.add_argument("--vmax", type=float, default=10, help="Maximum order of Bessel functions considered.")
parser.add_argument("--nOrder", type=int, default=20, help="Number of discrete orders of Bessel functions to consider.")
parser.add_argument("--zimin", type=float, default=-4,help="Minimum value of Im(z).")
parser.add_argument("--zimax", type=float, default=4, help="Maximum value of Im(z).")
parser.add_argument("--nzi", type=int, default=10, help="Number of points in Im(z)")
parser.add_argument("--zrmin", type=float, default=-10, help="Minimum value of Re(z).")
parser.add_argument("--zrmax", type=float, default=10, help="Maximum value of Re(z).")
parser.add_argument("--nzr", type=int, default=10, help="Number of points in  Re(z).")

args = parser.parse_args()

# -- MPI variables.
size = MPI.COMM_WORLD.size
rank = MPI.COMM_WORLD.rank 

# -- We verify that we have the proper number of processors.
if (args.pmax+1 != size):
  exit("The maximum order of differentiation must be equal to the number of processors.")

# ----------------- Data Declaration ------------------ #

# -- We open the file and create the datasets.
f = h5py.File(args.filename, 'w', driver='mpio', comm=MPI.COMM_WORLD)

datasets = [f.create_dataset("besselJ", (args.pmax+1, 2*args.nOrder+1,args.nzr,args.nzi),  dtype=complex),
            f.create_dataset("besselY", (args.pmax+1, 2*args.nOrder+1,args.nzr,args.nzi), dtype=complex),
            f.create_dataset("besselI", (args.pmax+1, 2*args.nOrder+1,args.nzr,args.nzi), dtype=complex),
            f.create_dataset("besselK", (args.pmax+1, 2*args.nOrder+1,args.nzr,args.nzi), dtype=complex),
            f.create_dataset("hankelH1", (args.pmax+1, 2*args.nOrder+1,args.nzr,args.nzi), dtype=complex),
            f.create_dataset("hankelH2", (args.pmax+1, 2*args.nOrder+1,args.nzr,args.nzi), dtype=complex)]

# -- We store the function callables in an array.
functions = [jvp, yvp, ivp, kvp, h1vp, h2vp]

# -- We prepare the values of the order and arguments.
orders = np.linspace(-args.vmax, args.vmax, 2*args.nOrder+1)
zr = np.linspace(args.zrmin, args.zrmax, args.nzr)
zi = np.linspace(args.zimin, args.zimax, args.nzi)


for i in range(len(datasets)):
  for j in range(2*args.nOrder+1):
    for k in range(args.nzr):
      for l in range(args.nzi):
        datasets[i].attrs.create("vmax", args.vmax)
        datasets[i].attrs.create("zimin", args.zimin)
        datasets[i].attrs.create("zimax", args.zimax)
        datasets[i].attrs.create("zrmin", args.zrmin)
        datasets[i].attrs.create("zrmax", args.zrmax)
        datasets[i][rank,j,k,l] = functions[i](orders[j], zr[k]+1j*zi[l], rank)

f.close()
