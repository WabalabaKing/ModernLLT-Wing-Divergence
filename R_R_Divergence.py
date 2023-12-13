# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 01:54:27 2023

@author: zheng
"""

from sympy import *
from scipy import optimize
import numpy as np
import LLTSolver 


####################################Aerodynamic definition
N = 10     #number of vorteex
sweep = -50  #LE sweep degree
b = 30       # halfspan 
Croot = 5       # chord
taper = 1  # taper ratio
Dihedral = 0  #degree
Cla = 2*np.pi*np.cos(sweep*np.pi/180)
AoA = 2
lam = sweep*np.pi/180
#Run to get CLa
cpy,Cl,CL,CD,Ar = LLTSolver.LLTSolver(N,sweep,b,Croot,taper,Dihedral,Cla,AoA)
CLa = Cl/(AoA*np.pi/180)
CLa = CLa[int(len(CLa)/2)::]
cpy = cpy[int(len(cpy)/2+0.5)::]
cpy=cpy-(cpy[1]-cpy[0])
Claa = np.polyfit(cpy, CLa, 3)
Ai = (cpy[1]-cpy[0])*Croot/np.cos(lam)
#################################Elastic Definition
GJ = 13330
EI = 8830
Order = 7
ECPy = np.linspace(0,b,Order)
e=0.25*Croot
#################################
x = symbols('y')
Clae = x**3*Claa[0]+x**2*Claa[1]+x*Claa[2]+Claa[3]
qd = symbols('qd',positive=True)
C = symbols('c0:%d'%(Order*2))
LHS = np.zeros((Order*2,Order*2))
RHS = np.zeros((Order*2,Order*2))
q1 = 0
q2=0
for i in range(Order):
    q1 = q1+C[i]*(x/b)**(i+2)   #theta
    q2 = q2+C[i+Order]*(x/b)**(i+2)   #dwdy
q1d = diff(q1,x) 
q2d = diff(q2,x) 
Utwist = 0.5*integrate(GJ*q1d**2,(x,0,b))
Ubend = 0.5*integrate(EI*q2d**2,(x,0,b))
Utot = Utwist+Ubend
dL =Clae*Ai/Ar*(AoA/180*np.pi+q1*np.cos(lam)-q2*np.cos(lam)*np.sin(lam))
Wtot = integrate(dL*q2,(x,0,b))+integrate(dL*e*q1,(x,0,b))
for i in range(Order*2):
   # RHS[i][i] = diff(Wtot,C[i])
    for j in range(Order*2):
        LHS[i][j] = diff(Utot,C[i],C[j])
        RHS[i][j] = diff(Wtot,C[i],C[j])
def addq(QD):
    Tot = LHS-RHS*QD
    Td = np.linalg.det(Tot)
    return Td

RigidQD = -EI*GJ/(GJ*np.cos(lam)**2*np.sin(lam)**2*b-EI*e*np.cos(lam)**3)
QD = optimize.fsolve(addq,10)
print(QD/12)
print(RigidQD/12)






# for i in range(Order):
#     Y = ECPy[i]
#     for j in range(Order):
#         LHS1[i][j] = EI*diff(q1ddd,C[j]).subs(x,Y)
#         LHS2[i+Order][j+Order] = GJ*diff(q2dd,C[j+Order]).subs(x,Y)
# LHS = LHS1+LHS2
# RHS1 = zeros(Order*2,Order*2)
# RHS2 = zeros(Order*2,Order*2)
# t = qd*Croot*e*CLa*(AoA*np.pi/180+q2*np.cos(lam)-q1*np.cos(lam)*np.sin(lam))
# dtdy = diff(t,x)
# z = qd*Croot*CLa*(AoA*np.pi/180+q2*np.cos(lam)-q1*np.cos(lam)*np.sin(lam))
# temp1 = z*np.cos(lam)+dtdy*np.cos(lam)*np.sin(lam)
# temp2 = t*np.cos(lam)*np.cos(lam)
# for i in range(Order):
#     Y = ECPy[i]
#     for j in range(Order*2):
#         RHS1[i,j] = diff(temp1,C[j]).subs(x,Y)
#         RHS2[i+Order,j] = diff(temp2,C[j]).subs(x,Y)
# RHS = RHS1+RHS2

# TLHS = RHS-LHS
# RigidQD = -EI*GJ/(GJ*np.cos(lam)**2*np.sin(lam)**2*b-EI*e*np.cos(lam)**3)
# QD = solve(TLHS.det(),qd)






