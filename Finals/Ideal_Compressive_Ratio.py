#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 09:56 2024

@author: nathanrawiri

The purpose of this file is to find a distribution of the strongest plys for compressive loads

"""



import math
import numpy as np
import matplotlib.pyplot as plt

import Build_ABD_Matrix as abdFile


import Material_Constants as cons
E1, E2, G12, v12, t, Xt, Xc, Yt, Yc, S, v21, Q_basic = cons.main() #Get material values


"""Initialise loads"""

Ny = 0; Ns = 0; Mx = 0; My = 0; Ms = 0;
Mz = 5 * 10**6 # MN

R = 3000

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

Ef1 = 120 #maybe change this
vf12 = 0.2  
mof = 1.1
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

#def getBucklingLoad(D, m, n, a, b, k): #b is width, a is length
#    D11 = D[0][0]; D12 = D[0][1]; D66 = D[2][2]; D22 = D[1][1]; AR = a/b
#    return np.pi * (D11*m**4 + 2 * (D12 + 2*D66)*m**2*n**2*AR**2+D22*n**4*AR**4)/(a**2*(m**2+k*n**2*AR**2))

def getBucklingLoadCCSS(D, m, a, b): #b is width, a is length
    D11 = D[0][0]; D12 = D[0][1]; D66 = D[2][2]; D22 = D[1][1]; lamb = (a/b)*(D22/D11)**(1/4)
    if (0 < lamb < 1.662):
        K = (4/lamb**2)+((2*(D12+2*D66))/(np.sqrt(D11*D22)))+((3/4)*lamb**2)
    elif (lamb >= 1.662):
        K = ((m**4+8*m**2+1)/(lamb**2*(m**2+1)))+((2*(D12+2*D66))/(np.sqrt(D11*D22)))+((lamb**2)/(m**2+1))
    return (np.pi**2/(b**2))*np.sqrt(D11*D22)*K
    


def getRandomLayup():
    return np.random.choice([0, 15, -15, 30,-30,45,-45,60,-60,75,-75,90,-90], size=(1, 6))

def getRandomPly():
    return np.random.choice([0, 15, -15, 30,-30,45,-45,60,-60,75,-75,90,-90])

def main():
    safe_ns = []
    base_layups = []
    bucklingStresses = []
    for i in range(0, 10000, 1):
        base_layup = getRandomLayup()[0]
        base_layup = base_layup.tolist()
        base_layups.append(base_layup)
        #layup = [getRandomPly()]
        """Start loop for getting the n at which it doesn't buckle"""
        this_failed = True;
        n_ = 1
        while this_failed:
            """Make ABD Matrix"""
            D, ABD, h = abdFile.getABDMatrix(base_layup, n_)                                                      #mm
    
            a = 500
            
            # b = 300 is good for this file because it cuts off the sample size 
            # to only the very best. 
            # If you change b to 450, you get a different distribution but
            # it's only because the sample size not tkaing only the very best
            # with 300, only the very best layups are selected.
            b = 300 
            m = 1
            
            Iz = np.pi*(R**4 - (R-h)**4)/4
            actualStress = Mz * R / Iz  # N/(mm^2)

            
            bucklingStress = getBucklingLoadCCSS(D, m, a, b) / (h)  # N/mm / mm = N / (mm^2)
            
            if bucklingStress > (2 * (actualStress)): #2x for factor of safety
                this_failed = False;

                safe_ns.append(n_)
                bucklingStresses.append(bucklingStress)
            n_ += 1

            

    min_n = np.min(safe_ns)
    min_weight_layups_w_load = [(base_layup, buckStress) for base_layup, n__, buckStress in zip(base_layups, safe_ns, bucklingStresses) if n__ == min_n]
    
    #print(values_at_indexes)
    ordered_lightest_layups = sorted(min_weight_layups_w_load, key=lambda x: x[1], reverse = True)
    #print(ordered_lightest_layups)
    print(f"Min N = {min_n}")
    print(f"num mins = {len(ordered_lightest_layups)}")

    value_counts = {}

    
    # Iterate over each list in values_at_indexes
    for sublist in ordered_lightest_layups:
        # Iterate over each value in the sublist
        for value in sublist[0]:
            # Check if the value is already in the dictionary
            if value in value_counts:
                # If it is, increment its count by 1
                value_counts[value] += 1
            else:
                # If it's not, initialize its count to 1
                value_counts[value] = 1
    
    # Extract values and counts from the dictionary
    values = list(value_counts.keys())
    counts = list(value_counts.values())
    
    # Create a dot plot
    plt.figure(figsize=(8, 6))
    plt.scatter(values, counts, color='blue')
    plt.xlabel('Angle of ply')
    plt.ylabel('Count')
    plt.title('Freq of angle present in lowest weight layups for compressive buckling')
    plt.grid(True)
    plt.xticks(values)
    plt.show()
    
    
    return 


main() 










