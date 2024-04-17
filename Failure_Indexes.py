#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 22:14:33 2024

@author: nathanrawiri
"""

import math
import numpy as np

#import Build_ABD_Matrix as abdFile
import Material_Constants as cons
E1, E2, G12, v12, t, Xt, Xc, Yt, Yc, S, v21, Q_basic = cons.main() #Get material values

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
    
"""For the strain per lamina matrix"""
def getStrainConversionMatrix(delta):
    array = [[math.cos(delta)**2, math.sin(delta)**2, math.sin(delta)*math.cos(delta)],
             [math.sin(delta)**2, math.cos(delta)**2, -math.sin(delta)*math.cos(delta)],
             [-2*math.sin(delta)*math.cos(delta), 2*math.sin(delta)*math.cos(delta), (math.cos(delta)**2-math.sin(delta)**2)]]
    return array
    

def checkFI(Nx_, Ns_, ABD, layup):   # It is ok to pass is the half symmetric layup to this because z is not used in the formula.
    angles = layup
    """Initialise loads"""
    Ny = 0; Mx = 0; My = 0; Ms = 0;
    Nx = Nx_; Ns = Ns_
    """Apply loads"""
    N = [[Nx],[Ny],[Ns]]
    M = [[Mx],[My],[Ms]]
    Nm = np.vstack((N, M))
    
    """ABD Inv"""
    ABD_inv = np.linalg.inv(ABD)
    
    """Get strain (e)"""
    ek = np.dot(ABD_inv, Nm) #strain and k matrix
    ex0 = ek[0][0]; ey0 = ek[1][0]; es0 = ek[2][0]  
    e0 = np.array([[ex0],[ey0],[es0]]) #mid-plain strain
    
    """Convert into lamina stress and strain"""
    strain = {}
    stress = {}
    FI = {}
    FI_max = 0
    puck_max = 0
    for i, angle in enumerate(angles):
        if angle is not None:
            """get stress and strain for angle"""
            strain[angle] = np.dot(getStrainConversionMatrix(np.deg2rad(angle)), e0)
            stress[angle] = np.dot(Q_basic, strain[angle])
            
            FI[angle] = [0,0,0]
            """Lamina directional stresses"""
            o_x = stress[angle][0][0]
            o_y = stress[angle][1][0]
            t_xy = stress[angle][2][0]
            """get max. failure index"""
            FI[angle][0] = o_x / Xt if o_x > 0 else o_x / -Xc
            FI[angle][1] = o_y / Yt if o_y > 0 else o_y / -Yc
            FI[angle][2] = abs(t_xy / S)

            Max = np.max(FI[angle])
            
            if Max > FI_max:
                FI_max = Max;
            
            """PUCK"""
            puckFF = getPuckFF(o_x, o_y)
            puckFI = getPuckFI(t_xy, o_y)
            puck = np.max([puckFF, puckFI])
            if puck > puck_max:
                puck_max = puck
            #print()
            
    FI_max = np.max((puck_max, FI_max))        
    return FI_max

"""
layup = [0, 0, -45, 45, 0, 0, 90]
D, ABD, h = abdFile.getABDMatrix(layup, 1)

FI_max = checkFI(531.0180783291642, 0, ABD, layup)
print(f"FI Max = {FI_max}")"""


