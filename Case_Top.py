#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 7 11:56 2024

@authors: Nathan Rawiri, Roman Crimi, James Bruce

Composites Assignement 1, Question 2a

[0/90/+-45]2s

Biaxial Stress Failure envelope for Puck and Max. failure criteria

Ny-Ns loading

"""



import math
import numpy as np
#import Ass2_Working_Puck as puck
import Build_ABD_Matrix as abd
import heapq
import matplotlib.pyplot as plt


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


#set_angles = [45, -45, 0, 0, 45, -45, 0, 0, 45, -45, 0, 0, 30,\
              #-30, 30, 0, 0, -30, 60, -60, 0, 0, 0, 90, 90, 0, 0, 90,  0, 0]
#set_angles = set_angles + set_angles[::-1]
#angles = set_angles
#h = len(angles) * t

Mz = 15e9                               # Nmm
R = 3000                                # mm
"""Iz = np.pi*(R**4 - (R-h)**4)/4          # mm^4

# LOAD  
bending_stress = Mz * R * t / Iz        # N / mm^2
arb_plate_width = 300                   #mm
Nx = bending_stress * arb_plate_width   # N / mm
print(f"Nx = {Nx}")"""

#other load
Ny = 0; Ns = 0; Mx = 0; My = 0; Ms = 0 


#thicknesses = [t] * len(angles)

v21 = v12 * E2/E1  #0.214
Q = 1-v12*v21
Q11 = E1 / Q 
Q22 = E2 / Q
Q12 = v12*E2 / Q
Q66 = G12 
Q_basic = [[Q11, Q12, 0], [Q12, Q22, 0],[0, 0, Q66]]


"""
Define Functions
"""

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

"""For the strain per lamina matrix"""
def getStrainConversionMatrix(delta):
    array = [[math.cos(delta)**2, math.sin(delta)**2, math.sin(delta)*math.cos(delta)],
             [math.sin(delta)**2, math.cos(delta)**2, -math.sin(delta)*math.cos(delta)],
             [-2*math.sin(delta)*math.cos(delta), 2*math.sin(delta)*math.cos(delta), (math.cos(delta)**2-math.sin(delta)**2)]]
    return array

def getRandomLayup():
    return np.random.choice([0, 15, -15, 30,-30,45,-45,60,-60,75,-75,90], size=(1, 6))

        

"""Initialise load, clean laminate etc. for each loading ratio"""

FI_maxs = []
layups = []
n = 1
for i in range(1, 50000, 1):
    layup = getRandomLayup()[0]
    ABD, h = abd.getABDMatrix(layup, n)
    ABD_inv = np.linalg.inv(ABD)
    layups.append(layup)
    
    Iz = np.pi*(R**4 - (R-h)**4)/4          # mm^4

    # LOAD  
    Nx = Mz * R * h / Iz        # N / mm                 #mm
    #print(f"Nx = {Nx}")
    
    
    
    
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
            
            FI[angle] = [0,0,0]
            """get failure index"""
            o_x = stress[angle][0][0]
            o_y = stress[angle][1][0]
            t_xy = stress[angle][2][0]
            FI[angle][0] = o_x / Xt if o_x > 0 else o_x / -Xc
            FI[angle][1] = o_y / Yt if o_y > 0 else o_y / -Yc
            FI[angle][2] = abs(t_xy / S)
            
            
            """ check if max failure index is 1"""
            Max = np.max(FI[angle])
            #if Max > 0.999:
            #    print("failed")
            if Max > FI_max:
                FI_max = Max;
    FI_maxs.append(FI_max)


"""Takes the layups with lowest FI's - Note, the more iterations, the better the sample"""
num_samples = 300 #change this to change the number of lowest layups
 
index_of_smallest_FIs = heapq.nsmallest(num_samples, enumerate(FI_maxs), key=lambda x: x[1])
index_of_smallest_FIs = sorted(index_of_smallest_FIs, key=lambda x: x[1], reverse = True)
#print(index_of_smallest_FIs)
layups_for_smallest_FIs_under_1 = [(layups[index], FI) for index, FI in index_of_smallest_FIs if FI < 1]
#print(layups_for_smallest_FIs_under_1)

"""Get the count of each ply in the list of the lowest FI layups"""
ply_counts = {}
# Iterate over each list in layups_for_smallest_FIs_under_1
for layup_and_FI in layups_for_smallest_FIs_under_1:
    # Iterate over each ply in the layup_and_FI
    print(f"FI_Max = {layup_and_FI[1]:.3f}, Layup = {layup_and_FI[0]}")
    for ply in layup_and_FI[0]:
        # Check if the ply is already in the dictionary
        if ply in ply_counts:
            # If it is, increment its count by 1
            ply_counts[ply] += 1
        else:
            # If it's not, initialize its count to 1
            ply_counts[ply] = 1

# Extract plys and counts from the dictionary
plys = list(ply_counts.keys())
counts = list(ply_counts.values())

# Create a dot plot
plt.figure(figsize=(8, 6))
plt.scatter(plys, counts, color='blue')
plt.xlabel('Angle of ply')
plt.ylabel('Count')
plt.title('Freq of angle present in lowest FI layups for tensile Max')
plt.grid(True)
plt.xticks(plys)
plt.ylim(ymin=0)#, ymax=150)  # Adjust the limits for the y axis
plt.show()









