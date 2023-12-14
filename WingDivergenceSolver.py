# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 00:07:23 2023

@author: zheng
"""

from sympy import *
from scipy import optimize
import numpy as np
import LLTSolver 
def WingDivergenceSolver(N,sweep,b,Croot,taper,Dihedral,Cla,AoA,Order,GJ,EI,e):
    lam = sweep*np.pi/180
    #Run to get CLa
    cpy,Cl,CL,CD,Ar = LLTSolver.LLTSolver(N,sweep,b,Croot,taper,Dihedral,Cla,AoA)
    CLa = Cl/(AoA*np.pi/180)
    CLa = CLa[int(len(CLa)/2)::]
    cpy = cpy[int(len(cpy)/2+0.5)::]
    cpy=cpy-(cpy[1]-cpy[0])
    Claa = np.polyfit(cpy, CLa, 3)
    Ai = (cpy[1]-cpy[0])*Croot/np.cos(lam)
    ECPy = np.linspace(0,b,Order)
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
        q1 = q1+C[i]*(x/b)**(i+3)   #theta
        q2 = q2+C[i+Order]*(x/b)**(i+3) #dw
    q1d = diff(q1,x) 
    q2d = diff(q2,x) 
    q2dd = diff(q2,x,x) 
    Utwist = 0.5*integrate(GJ*q1d**2,(x,0,b))
    Ubend = 0.5*integrate(EI*q2dd**2,(x,0,b))
    Utot = Utwist+Ubend
    dL =Clae*Ai/Ar*(AoA/180*np.pi*np.cos(lam)+q1*np.cos(lam)-q2d*np.cos(lam)*np.sin(lam))
    Wtot = integrate(dL*q2d,(x,0,b))+integrate(dL*e*q1d,(x,0,b))
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
    print(QD)
    return QD