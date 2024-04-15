#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:13:37 2024

@author: nathanrawiri
"""

import math
import numpy as np

def m_(deg):
    return math.cos(np.deg2rad(deg))

def n_(deg):
    return math.sin(np.deg2rad(deg))


def Qxx_(d):
    m = m_(d); n = n_(d);
    return Q11 * m**4 + 2* (Q12 + 2*Q66) * m**2 * n**2 + Q22 * n**4

def Qxy_(d):
    m = m_(d); n = n_(d);
    return (Q11 + Q22 - 4*Q66)*m**2*n**2 + Q12*(m**4+n**4)

def Qyy_(d):
    m = m_(d); n = n_(d);
    return Q11*n**4 + 2*(Q12+2*Q66)*m**2*n**2 + Q22*m**4

def Qxs_(d):
    m = m_(d); n = n_(d);
    return ((Q11-Q12-2*Q66)*n*m**3 + (Q12-Q22+2*Q66)*n**3*m)

def Qys_(d):
    m = m_(d); n = n_(d);
    return ((Q11-Q12-2*Q66)*m*n**3 + (Q12-Q22+2*Q66)*m**3*n)

def Qss_(d):
    m = m_(d); n = n_(d);
    return (Q11+Q22-2*Q12-2*Q66)*n**2*m**2 + Q66*(n**4+m**4)



"""Assignment Values"""
E1 = 140
E2 = 11.2  
G12 = 5
v12 = 0.33

t = 0.135
#Need strength
Xt = 2200 # MPa
Xc = 1800 #MPa
Yt = 78 #MPa
Yc = 300 #MPa
S = 100 #MPa


v21 = v12 * E2/E1  #0.214
Q = 1-v12*v21
Q11 = E1 / Q 
Q22 = E2 / Q
Q12 = v12*E2 / Q
#Q21 =#does this equal Q12 every every time???
Q66 = G12
    
def getABDMatrix(angles, symmetry):
    
    #angles = [15, delta, -delta, 75, 75]
    angles = ( angles + angles[::-1] ) * symmetry 
    #print(angles)
    
    thicknesses = [0.135] * len(angles)
    
    """convert thickness into z values"""
    h = sum(thicknesses)
    z = []
    z.append(-h/2)
    h_i = 0
    for t in thicknesses:
        h_i += t
        z.append(h_i - (h/2)) #Maybe add abs to this if needed?
    #print(f"z = {z}")
    
    
    Q_matrix = {}
    
    """Create Q matrix for each angle (currenlty repeats if repetitive angle)"""
    for angle in angles:
        Qxx = Qxx_(angle)
        Qxy = Qxy_(angle)
        Qyy = Qyy_(angle)
        Qxs = Qxs_(angle)
        Qys = Qys_(angle)
        Qss = Qss_(angle)
        
        Q_matrix[angle] = [[Qxx, Qxy, Qxs], [Qxy, Qyy, Qys], [Qxs, Qys, Qss]]
    
    """Make ABD Matrix"""
    A = np.zeros((3,3))
    B = np.zeros((3,3))
    D = np.zeros((3,3))
    for k, angle in enumerate(angles):
        if angle is not None:
            for i in range(3):
                for j in range(3):
                    A[i][j] += Q_matrix[angle][i][j] * (z[k+1] - z[k])
                    B[i][j] += 1/2 * Q_matrix[angle][i][j] * (z[k+1]**2 - z[k]**2)
                    D[i][j] += 1/3 * Q_matrix[angle][i][j] * (z[k+1]**3 - z[k]**3)
    AB = np.hstack((A, B))
    BD = np.hstack((B, D))
    ABD = np.vstack((AB, BD))
    
    return ABD, h