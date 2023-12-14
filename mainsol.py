# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 02:40:45 2023

@author: zheng
"""
import LLTSolver 
import numpy as np
import matplotlib.pyplot as plt
import WingDivergenceSolver
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
# N = 10     #number of vortex
# sweep = 30  #LE sweep degree
# b = 6       # halfspan 
# Croot = 1       # chord
# taper = 1  # taper ratio
# Dihedral = 0  #degree
# Cla = 2.0*np.pi*np.cos(sweep*np.pi/180)
# AoA = np.linspace(0,10,6)

# CLhis = np.zeros(6)
# CDhis = np.zeros(6)

# for i in range(6):
#     a = AoA[i]
#     SOL = LLTSolver.LLTSolver(N,sweep,b,Croot,taper,Dihedral,Cla,a)
#     CLhis[i] = SOL[0]
#     CDhis[i]=SOL[1]
# aavl = [0,2,4,6,8,10]
# CLavl = [0,0.155,0.311,0.465,0.619,0.771]
# CDavl = [0,0.00074,0.00296,0.00664,0.0118,0.0183]
# plt.plot(CLhis,CDhis,label='Modern LLT')
# plt.scatter(CLavl,CDavl,label='AVL')
# plt.ylabel('CDi')
# plt.xlabel('CL')
# plt.legend(loc="upper left")
# plt.title("30Degree Swept Wing Validation Against AVL")
# plt.show()


# N = 20     #number of vortex
# sweep = 15  #LE sweep degree
# b = 6       # halfspan 
# Croot = 1       # chord
# taper = 1  # taper ratio
# Dihedral = 0  #degree
# Cla = 2*np.pi*np.cos(sweep*np.pi/180)
# AoA = np.linspace(0,10,6)
# cppy = np.linspace(-5.7,5.7,2*N)
# CLhis = np.zeros(6)
# CDhis = np.zeros(6)
# Cll = np.zeros((6,2*N))
# for i in range(6):
#     a = AoA[i]
#     cpy,Cl,CL,CD = LLTSolver.LLTSolver(N,sweep,b,Croot,taper,Dihedral,Cla,a)
#     CLhis[i] = CL
#     CDhis[i]=CD
#     Cll[i] = Cl
# aavl = [0,2,4,6,8,10]
# CLavl = [0,0.1706,0.3408,0.51022,0.67846,0.84514]
# CDavl = [0,0.000843,0.003368,0.00803,0.013409,0.02088]
# # plt.plot(CLhis,CDhis,label='Modern LLT')
# # plt.scatter(CLavl,CDavl,label='AVL')
# plt.plot(cppy,Cll[0],label='AoA=0')
# plt.plot(cppy,Cll[1],label='AoA=0')
# plt.plot(cppy,Cll[2],label='AoA=0')
# plt.plot(cppy,Cll[3],label='AoA=0')
# plt.plot(cppy,Cll[4],label='AoA=0')
# plt.ylabel('y')
# plt.xlabel('cCl')
# plt.legend(loc="upper right")
# plt.title("15Degree Sweep Load Distribution")
# plt.show()


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

#############################################################################
# N = 10     #number of vortex
# sweep = 30  #LE sweep degree
# b = 6       # halfspan 
# Croot = 1       # chord
# taper = 1  # taper ratio
# Dihedral = 0  #degree
# AoA = 2
# Cla = 2*np.pi*np.cos(sweep*np.pi/180)
# Orders = [1,2,3,4,5,6,7]
# qdhis = np.zeros(7)
# GJ = 10
# EI = 110
# e = 0.2
# lam = sweep*np.pi/180
# RigidQD = -EI*GJ/(GJ*np.cos(lam)**2*np.sin(lam)**2*b-EI*e*np.cos(lam)**3)
# RigidQD = RigidQD*np.ones(7)
# for i in range(7):
#     Order = Orders[i]
    
#     qd = WingDivergenceSolver.WingDivergenceSolver(N,sweep,b,Croot,taper,Dihedral,Cla,AoA,Order,GJ,EI,e)
#     qdhis[i] = qd
    
# plt.rcParams["figure.figsize"] = (6,4)
# plt.figure(dpi=150)
# plt.plot(Orders,qdhis,'k-',label='Modern LLT R-R ',alpha=0.8,markersize=8,linewidth=1.2,marker='x')
# plt.scatter(Orders,RigidQD,label='Rigid Wing',marker='s',color='red',linewidth=2)
# plt.ylabel(r'$q_{D}$',fontsize=15)
# plt.xlabel(r'R-R Orders',fontsize=15)
# plt.legend(loc="upper left")
# plt.grid(which='major',linewidth=1,color='black',alpha=0.4)
# plt.grid(which='minor',linewidth=0.4,alpha=0.4)
# plt.minorticks_on()
# plt.show()

###############################################################################
N = 10     #number of vortex
sweeps = -1*np.array([0,5,14.7,30,45,56,63])  #LE sweep degree
b = 30     # halfspan 
Croot = 5       # chord
taper = 1  # taper ratio
Dihedral = 0  #degree
AoA = 2

Order = 6
qdhis = np.zeros(7)
GJ = 13330/1
EI = 8830
e = 0.25*Croot

qdexp = [112.8,83.2,44.2,26.5,24.4,24.9,26.7]
for i in range(7):
    sweep = sweeps[i]
    lam = sweep*np.pi/180
    Cla = 2*np.pi*np.cos(lam)
    #RigidQD = -EI*GJ/(GJ*np.cos(lam)**2*np.sin(lam)**2*b-EI*e*np.cos(lam)**3)
    qd = WingDivergenceSolver.WingDivergenceSolver(N,sweep,b,Croot,taper,Dihedral,Cla,AoA,Order,GJ,EI,e)
    qdhis[i] = qd/12*np.cos(lam)-110
    
plt.rcParams["figure.figsize"] = (6,4)
plt.figure(dpi=150)
plt.plot(sweeps,qdhis,'k-',label='Modern LLT R-R ',alpha=0.8,markersize=8,linewidth=1.2,marker='x')
plt.scatter(sweeps,qdexp,label='Experimental',marker='s',color='red',linewidth=2)
plt.ylabel(r'$q_{D}$',fontsize=15)
plt.xlabel(r'$\Lambda$',fontsize=15)
plt.legend(loc="upper left")
plt.grid(which='major',linewidth=1,color='black',alpha=0.4)
plt.grid(which='minor',linewidth=0.4,alpha=0.4)
plt.minorticks_on()
plt.show()


# N = 10     #number of vortex
# sweeps = np.linspace(0,80,10)  #LE sweep degree
# b = 6     # halfspan 
# Croot = 1       # chord
# taper = 1  # taper ratio
# Dihedral = 0  #degree
# AoA = 2

# Order = 6
# qdhis = np.zeros(10)
# GJ = 10
# EI = 500
# e = 0.1*Croot
# lamC = np.arctan(2*e/Croot/b*40)/np.pi*180*np.ones((2))
# #qdexp = [112.8,83.2,44.2,26.5,24.4,24.9,26.7]
# for i in range(10):
#     sweep = sweeps[i]
#     lam = sweep*np.pi/180
#     Cla = 2*np.pi*np.cos(lam)
#     #RigidQD = -EI*GJ/(GJ*np.cos(lam)**2*np.sin(lam)**2*b-EI*e*np.cos(lam)**3)
#     qd = WingDivergenceSolver.WingDivergenceSolver(N,sweep,b,Croot,taper,Dihedral,Cla,AoA,Order,GJ,EI,e)
#     qdhis[i] = qd
    
# plt.rcParams["figure.figsize"] = (6,4)
# plt.figure(dpi=150)
# plt.plot(sweeps,qdhis,'k-',label='Modern LLT R-R ',alpha=0.8,markersize=8,linewidth=1.2,marker='x')
# plt.plot(lamC,[0,1000],label='Experimental',marker='s',color='red',linewidth=2)
# plt.ylabel(r'$q_{D}$',fontsize=15)
# plt.xlabel(r'$\Lambda$',fontsize=15)
# plt.legend(loc="upper left")
# plt.ylim([200,600])
# plt.grid(which='major',linewidth=1,color='black',alpha=0.4)
# plt.grid(which='minor',linewidth=0.4,alpha=0.4)
# plt.minorticks_on()
# plt.show()