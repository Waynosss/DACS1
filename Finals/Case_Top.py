#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 7 11:56 2024

Description: This file attempts to find the distribution of best angle plys for the top case, in tension.

@authors: Nathan Rawiri, Roman Crimi, James Bruce


"""



import math
import numpy as np
#import Ass2_Working_Puck as puck
import Build_ABD_Matrix as abd
import heapq
import matplotlib.pyplot as plt

import Material_Constants as cons
E1, E2, G12, v12, t, Xt, Xc, Yt, Yc, S, v21, Q_basic = cons.main() #Get material values



Mz = 15e9                               # Nmm
R = 3000                                # mm

#other load
Ny = 0; Ns = 0; Mx = 0; My = 0; Ms = 0 




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
for i in range(0, 100000, 1):
    layup = getRandomLayup().tolist()[0]
    D, ABD, h = abd.getABDMatrix(layup, n)
    ABD_inv = np.linalg.inv(ABD)
    layups.append(layup)
    
    Iz = np.pi*(R**4 - (R-h)**4)/4          # mm^4

    # LOAD  
    Nx = Mz * R * h / Iz        # N / mm                 #mm
    #print(f"Nx = {Nx}")         # N / mm
    
    
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
    #print(FI_max)
    #print(layup)



"""Takes the layups with lowest FI's - Note, the more iterations, the better the sample"""
num_samples = 300 #change this to change the number of lowest layups
 
index_of_smallest_FIs = heapq.nsmallest(num_samples, enumerate(FI_maxs), key=lambda x: x[1])
index_of_smallest_FIs = sorted(index_of_smallest_FIs, key=lambda x: x[1], reverse = True)
#print(index_of_smallest_FIs)
layups_for_smallest_FIs_under_1 = [(layups[index], FI) for index, FI in index_of_smallest_FIs]# if FI < 1]
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









