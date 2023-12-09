# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 02:38:20 2023

@author: zheng
"""
import numpy as np
from scipy import optimize
def LLTSolver(N,sweep,b,Croot,taper,Dihedral,Cla,AoA):
    Ctip = Croot*taper
    sweep = sweep/360*2*np.pi  #deg2rad
    Dihedral = Dihedral/360*2*np.pi  #deg2rad
    Ar = (Croot+Ctip)*b
    #################Definition of Airfoil and Flow

    AoA = AoA/360*2*np.pi

    ui =np.array( [np.cos(AoA),0,np.sin(AoA)])
    ui = ui/np.linalg.norm(ui)
    un = np.array( [-1*np.sin(AoA),0,np.cos(AoA)])
    un = un/np.linalg.norm(un)
    #####################Geometry discretization###########

    LETip = np.array([b*np.tan(sweep),b,b*np.tan(Dihedral)])
    cPTip = np.array([b*np.tan(sweep)+0.25*Ctip,b,b*np.tan(Dihedral)])
    LERoot = np.array([0,0,0])
    cPRoot = np.array([0.25*Croot,0,0])
    cPX = np.linspace(cPRoot[0],cPTip[0],N+1)
    cPY = np.linspace(cPRoot[1],cPTip[1],N+1)
    cPZ = np.linspace(cPRoot[2],cPTip[2],N+1)
    cib = np.linspace(Croot,Ctip,N)

    cPX = np.concatenate((np.flip(cPX[1::]),cPX),axis=None)
    cPY = np.concatenate((-1*np.flip(cPY[1::]),cPY),axis=None)
    cPZ = np.concatenate((np.flip(cPZ[1::]),cPZ),axis=None)
    cib = np.concatenate((np.flip(cib),cib),axis=None)

    #######################LHS loop
    def LLT(G):
        coeff = np.empty(2*N)
        
        for j in range(2*N):
            VJ = np.zeros([1,3])
            CPj = np.array([(cPX[j]+cPX[j+1])/2,(cPY[j]+cPY[j+1])/2,(cPZ[j]+cPZ[j+1])/2])
            r0 = np.array([cPX[j+1],cPY[j+1],cPZ[j+1]])-np.array([cPX[j],cPY[j],cPZ[j]]);
            dl = r0/np.linalg.norm(r0)
            for i in range(2*N):
               
                ri2j = CPj- np.array([cPX[i+1],cPY[i+1],cPZ[i+1]])
                ri1j = CPj- np.array([cPX[i],cPY[i],cPZ[i]])

                
                ri2jn = np.linalg.norm(ri2j)
                ri1jn = np.linalg.norm(ri1j)
                if i==j:
                    temp= cib[i]/(4*np.pi) *(\
                             np.cross(ui,ri2j)/(ri2jn*(ri2jn-np.dot(ui,ri2j)))- \
                             np.cross(ui,ri1j)/(ri1jn*(ri1jn-np.dot(ui,ri1j))) )  
                    
                else:
                    temp = cib[i]/(4*np.pi) *(\
                             ((ri1jn+ri2jn)*np.cross(ri1j,ri2j))/(ri1jn*ri2jn*(ri1jn*ri2jn+np.dot(ri1j,ri2j)))+\
                             np.cross(ui,ri2j)/(ri2jn*(ri2jn-np.dot(ui,ri2j)))- \
                             np.cross(ui,ri1j)/(ri1jn*(ri1jn-np.dot(ui,ri1j))) )  
                VJ=VJ+temp*G[i]
            coeffs =2*np.cross((VJ+ui),dl)
            coeffL = 2*np.linalg.norm(coeffs)*G[j]
           # rhs = Cla*np.dot(np.arctan(ui+VJ),np.array([np.sin(AoA),0,np.cos(AoA)]))\
           #     /np.dot(np.arctan(ui+VJ),np.array([np.cos(AoA),0,-1*np.sin(AoA)]))
            ai = np.arctan2(np.dot(ui+VJ,np.array([np.sin(AoA),0,np.cos(AoA)])),np.dot(ui+VJ,np.array([np.cos(AoA),0,-1*np.sin(AoA)])))
            coeff[j] = coeffL-ai*Cla
        return coeff


    solG = optimize.root(LLT,np.ones(2*N)*0.01)
    Gamma = solG.x          

    #Solve for CL
    CF=np.zeros(3)
    for i in range(2*N):
        CPi = np.array([(cPX[i]+cPX[i+1])/2,(cPY[i]+cPY[i+1])/2,(cPZ[i]+cPZ[i+1])/2])
        r0 = np.array([cPX[i+1],cPY[i+1],cPZ[i+1]])-np.array([cPX[i],cPY[i],cPZ[i]]);
        dl = r0/np.linalg.norm(r0)
        temp=0
        for j in range(2*N):
            rj2i = CPi- np.array([cPX[j+1],cPY[j+1],cPZ[j+1]])
            rj1i = CPi- np.array([cPX[j],cPY[j],cPZ[j]])
            rj2in = np.linalg.norm(rj2i)
            rj1in = np.linalg.norm(rj1i)
            if i==j:
                vji= cib[i]/(4*np.pi) *(\
                         np.cross(ui,rj2i)/(rj2in*(rj2in-np.dot(ui,rj2i)))- \
                         np.cross(ui,rj1i)/(rj1in*(rj1in-np.dot(ui,rj1i))) )  
                
            else:
                vji = cib[i]/(4*np.pi) *(\
                         ((rj1in+rj2in)*np.cross(rj1i,rj2i))/(rj1in*rj2in*(rj1in*rj2in+np.dot(rj1i,rj2i)))+\
                         np.cross(ui,rj2i)/(rj2in*(rj2in-np.dot(ui,rj2i)))- \
                         np.cross(ui,rj1i)/(rj1in*(rj1in-np.dot(ui,rj1i))) )
            temp = temp+Gamma[i]*Gamma[j]*vji
        CF =CF+ np.cross((Gamma[i]*ui+temp),dl)*(cib[i]*(cPY[i+1]-cPY[i]))/Ar
    CF=CF*2
    CL = np.linalg.norm(np.dot(CF,un)*un)
    CDi = np.linalg.norm(np.dot(CF,ui)*ui)
    print(CL)
    return [CL,CDi]