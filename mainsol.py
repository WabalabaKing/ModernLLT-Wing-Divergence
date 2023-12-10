# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 02:40:45 2023

@author: zheng
"""
import LLTSolver 
import numpy as np
import matplotlib.pyplot as plt

#################################StraightWing Calibration########################
# N = 10     #number of vortex
# sweep = 0  #LE sweep degree
# b = 6       # halfspan 
# Croot = 1       # chord
# taper = 1  # taper ratio
# Dihedral = 0  #degree
# Cla = 1.9*np.pi*np.cos(sweep)
# AoA = np.linspace(0,10,6)

# CLhis = np.zeros(6)
# CDhis = np.zeros(6)

# for i in range(6):
#     a = AoA[i]
#     SOL = LLTSolver.LLTSolver(N,sweep,b,Croot,taper,Dihedral,Cla,a)
#     CLhis[i] = SOL[0]
#     CDhis[i]=SOL[1]
# aavl = [0,2,4,6,8,10]
# CLavl = [0,0.175,0.350,0.52408,0.69695,0.86827]
# CDavl = [0,0.0008598,0.0034,0.00767,0.01367,0.0213]
# plt.plot(CLhis,CDhis,label='Modern LLT')
# plt.scatter(CLavl,CDavl,label='AVL')
# plt.ylabel('CDi')
# plt.xlabel('CL')
# plt.legend(loc="upper left")
# plt.title("Straight Wing Validation Against AVL")
# plt.show()

#################################SweptWing Calibration########################
N = 10     #number of vortex
sweep = 30  #LE sweep degree
b = 6       # halfspan 
Croot = 1       # chord
taper = 1  # taper ratio
Dihedral = 0  #degree
Cla = 2.0*np.pi*np.cos(sweep*np.pi/180)
AoA = np.linspace(0,10,6)

CLhis = np.zeros(6)
CDhis = np.zeros(6)

for i in range(6):
    a = AoA[i]
    SOL = LLTSolver.LLTSolver(N,sweep,b,Croot,taper,Dihedral,Cla,a)
    CLhis[i] = SOL[0]
    CDhis[i]=SOL[1]
aavl = [0,2,4,6,8,10]
CLavl = [0,0.155,0.311,0.465,0.619,0.771]
CDavl = [0,0.00074,0.00296,0.00664,0.0118,0.0183]
plt.plot(CLhis,CDhis,label='Modern LLT')
plt.scatter(CLavl,CDavl,label='AVL')
plt.ylabel('CDi')
plt.xlabel('CL')
plt.legend(loc="upper left")
plt.title("30Degree Swept Wing Validation Against AVL")
plt.show()

# N = 5     #number of vortex
# sweep = 60  #LE sweep degree
# b = 6       # halfspan 
# Croot = 1       # chord
# taper = 1  # taper ratio
# Dihedral = 0  #degree
# Cla = 2*np.pi
# AoA = np.linspace(0,10,6)

# CLhis = np.zeros(6)
# CDhis = np.zeros(6)

# for i in range(6):
#     a = AoA[i]
#     SOL = LLTSolver.LLTSolver(N,sweep,b,Croot,taper,Dihedral,Cla,a)
#     CLhis[i] = SOL[0]
#     CDhis[i]=SOL[1]
# aavl = [0,2,4,6,8,10]
# CLavl = [0,0.096,0.192,0.289,0.384,0.478]
# CDavl = [0,0.00034,0.00134,0.00301,0.00533,0.0082966]
# plt.plot(CLhis,CDhis,label='Modern LLT')
# plt.scatter(CLavl,CDavl,label='AVL')
# plt.ylabel('CDi')
# plt.xlabel('CL')
# plt.legend(loc="upper right")
# plt.title("60Degree Swept Wing Validation Against AVL")
# plt.show()