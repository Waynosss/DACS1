#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 7 11:56 2024

Description: This file analyses the selected layup for the side of the fuselage

@authors: Nathan Rawiri, Roman Crimi, James Bruce


"""


import math
import numpy as np
#import Ass2_Working_Puck as puck
import Build_ABD_Matrix as abd



import Material_Constants as cons
E1, E2, G12, v12, t, Xt, Xc, Yt, Yc, S, v21, Q_basic = cons.main() #Get material values



Mz = 15e9                               # Nmm
V = 1.5e6                               # N
R = 3000                                # mm

#other load
Ny = 0; Mx = 0; My = 0; Ms = 0 



"""For the strain per lamina matrix"""
def getStrainConversionMatrix(delta):
    array = [[math.cos(delta)**2, math.sin(delta)**2, math.sin(delta)*math.cos(delta)],
             [math.sin(delta)**2, math.cos(delta)**2, -math.sin(delta)*math.cos(delta)],
             [-2*math.sin(delta)*math.cos(delta), 2*math.sin(delta)*math.cos(delta), (math.cos(delta)**2-math.sin(delta)**2)]]
    return array

"""Puck"""
p12p = 0.3
p12n = 0.2
o23A = (S / (2*p12n)) * (math.sqrt(1 + 2*p12n * (Yc/S))-1)
p23n = p12n * o23A / S
def modeA(t12, o2):
    return math.sqrt((t12 / S)**2 + (1 - p12p*(Yt / S))**2 * (o2 / Yt)**2 ) + p12p * o2 / S

def modeB(t12, o2):
    return (1/S) * (math.sqrt(t12**2 + (p12n*o2)**2) + p12n*o2)

def modeC(t12, o2):
    return ((t12 / (2*(1 + p23n) * S))**2 + (o2/Yc)**2)*(Yc/(-o2))

Ef1 = 225
vf12 = 0.2  
mof = 1.1
"""Define Puck Functions"""
def getPuckFF(o1, o2):
    R = Xc if o1 < 0 else Xt
    return (1/R) * (o1 - (v12-vf12*mof*E1/Ef1)*o2)
    
def getPuckFI(t12, o2):
    #print(f"Puck sigma2 = {o2}")
    o23_max = abs(o23A/abs(S))
    mid = abs(o2/t12)
    if o2 > 0:
        return modeA(t12, o2)
    elif (0 <= mid <= o23_max):
        return modeB(t12, o2)
    else:
        return modeC(t12, o2)


def getRandomLayup():
    return np.random.choice([0, 15, -15, 30,-30,45,-45,60,-60,75,-75,90], size=(1, 6))


"""Initialise load, clean laminate etc. for each loading ratio"""

layups = []
n = 1
#layup = [15, 0, 0, -15, -15, 0]#FI_Max = 0.050, Layup = [15, 0, 0, -15, -15, 0]
layup = [0, 45, -45, 90]
D, ABD, h = abd.getABDMatrix(layup, n)
ABD_inv = np.linalg.inv(ABD)
layups.append(layup)

d = np.linalg.inv(D)
E1b = 12 / (h**3 * d[0][0])
#print(f"E1b = {E1b}")

Iz = np.pi*(R**4 - (R-h)**4)/4          # mm^4

A = np.pi * (R**2 - (R-h)**2)
Tau_max = 2 * V / A

# LOAD  
#Nx = Mz * R * h / Iz        # N / mm                 #mm
Nx= 0
Ns = Tau_max * h
print(f"Nxy = {Ns}")

print(f'thickness = {h} mm')


"""Apply loads"""
N = [[Nx],[Ny],[Ns]]
M = [[Mx],[My],[Ms]]
Nm = np.vstack((N, M))

"""Get strain (e)"""
ek = np.dot(ABD_inv, Nm) #strain and k matrix
ex0 = ek[0][0]; ey0 = ek[1][0]; es0 = ek[2][0]  
e0 = np.array([[ex0],[ey0],[es0]]) #mid-plain strain

"""Get stress (o)"""
o0 = np.dot(Q_basic, e0);  
oy = o0[1]
t12 = o0[2]

"""Convert into lamina stress and strain"""
strain = {}
stress = {}
FI = {}
FI_max = 0
for i, angle in enumerate(layup):
    Max = 0
    if angle is not None:
        """get stress and strain for angle"""
        strain[angle] = np.dot(getStrainConversionMatrix(np.deg2rad(angle)), e0)
        stress[angle] = np.dot(Q_basic, strain[angle])
        
        FI[angle] = [0,0,0, 0, 0]
        """get failure index"""
        o_x = stress[angle][0][0]
        o_y = stress[angle][1][0]
        t_xy = stress[angle][2][0]
        FI[angle][0] = o_x / Xt if o_x > 0 else o_x / -Xc
        FI[angle][1] = o_y / Yt if o_y > 0 else o_y / -Yc
        FI[angle][2] = abs(t_xy / S)
        
        puckFI = getPuckFI(t_xy, o_y)
        puckFF = getPuckFF(o_x, o_y)
        
        FI[angle][3] = puckFF 
        FI[angle][4] = puckFI
        """ check if max failure index is 1"""
        Max = np.max(FI[angle])

        if Max > FI_max:
            FI_max = Max;
        
        
        
        
print(f"Failure Index = {FI_max}")

print(f"Layup : {layup}")

print("Note this is the minimum possible layup that follows requirements / standard practice. But it is enough to withstand the shear force")



