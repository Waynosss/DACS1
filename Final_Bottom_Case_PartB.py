#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 09:56 2024

@author: nathanrawiri


"""



import math
import numpy as np
import matplotlib.pyplot as plt
import heapq
import random
import Material_Constants as cons

E1, E2, G12, v12, t, Xt, Xc, Yt, Yc, S, v21, Q_basic = cons.main() #Get material values


#set_angles = [0, 90, 45, -45, -45, 45, 90, 0, 0, 90, 45, -45, -45, 45, 90, 0]
#angles = set_angles


"""Initialise loads"""

Ny = 0; Ns = 0; Mx = 0; My = 0; Ms = 0;
Mz = 15 * 10**6 # MN


R = 3000


def m_(deg):
    return math.cos(np.deg2rad(deg))

def n_(deg):
    return math.sin(np.deg2rad(deg))

"""Change all Q formulas to use a Q matrix input"""
def Qxx_(d, Q):
    m = m_(d); n = n_(d);
    return Q[0][0] * m**4 + 2* (Q[0][1] + 2*Q[2][2]) * m**2 * n**2 + Q[1][1] * n**4

def Qxy_(d, Q):
    m = m_(d); n = n_(d);
    return (Q[0][0] + Q[1][1] - 4*Q[2][2])*m**2*n**2 + Q[0][1]*(m**4+n**4)

def Qyy_(d, Q):
    m = m_(d); n = n_(d);
    return Q[0][0]*n**4 + 2*(Q[0][1]+2*Q[2][2])*m**2*n**2 + Q[1][1]*m**4

def Qxs_(d, Q):
    m = m_(d); n = n_(d);
    return ((Q[0][0]-Q[0][1]-2*Q[2][2])*n*m**3 + (Q[0][1]-Q[1][1]+2*Q[2][2])*n**3*m)

def Qys_(d, Q):
    m = m_(d); n = n_(d);
    return ((Q[0][0]-Q[0][1]-2*Q[2][2])*m*n**3 + (Q[0][1]-Q[1][1]+2*Q[2][2])*m**3*n)

def Qss_(d, Q):
    m = m_(d); n = n_(d);
    return (Q[0][0]+Q[1][1]-2*Q[0][1]-2*Q[2][2])*n**2*m**2 + Q[2][2]*(n**4+m**4)

"""For the strain per lamina matrix"""
def getStrainConversionMatrix(delta):
    array = [[math.cos(delta)**2, math.sin(delta)**2, math.sin(delta)*math.cos(delta)],
             [math.sin(delta)**2, math.cos(delta)**2, -math.sin(delta)*math.cos(delta)],
             [-2*math.sin(delta)*math.cos(delta), 2*math.sin(delta)*math.cos(delta), (math.cos(delta)**2-math.sin(delta)**2)]]
    return array


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

def getBucklingLoadSSSS(D, m, a, b):
    D11 = D[0][0]; D12 = D[0][1]; D66 = D[2][2]; D22 = D[1][1];
    #D11 = 659.7; D12 = 466.9; D22 = 659.7; D66 = 494;
    AR = a/b;
    return np.pi**2 * (D11*m**4 + 2*(D12+2*D66)*m**2*AR**2 + D22*AR**4 ) / (a**2 * m**2)
    
"""ABD Matrix Function"""
def getABDMatrix(layup, symmetry):
    
    #angles = [0, 90, 45, -45, delta]
    """The below line creates the layup with a symmetric laminate * n symmetry"""
    angles = layup
    angles = ( angles + angles[::-1] ) * symmetry    

    thicknesses = [0.135] * len(angles)
    """convert thickness into z values"""
    h = sum(thicknesses)
    z = []
    z.append(-h/2)
    h_i = 0
    for t in thicknesses:
        h_i += t
        z.append(h_i - (h/2))
    
    Q_matrix = {}
    """Create Q matrix for each angle (currenlty repeats if repetitive angle)"""
    for angle in angles:
        Qxx = Qxx_(angle, Q_basic)
        Qxy = Qxy_(angle, Q_basic)
        Qyy = Qyy_(angle, Q_basic)
        Qxs = Qxs_(angle, Q_basic)
        Qys = Qys_(angle, Q_basic)
        Qss = Qss_(angle, Q_basic)
        
        Q_matrix[angle] = [[Qxx, Qxy, Qxs], [Qxy, Qyy, Qys], [Qxs, Qys, Qss]]
    
    """Make ABD Matrix"""
    A = np.zeros((3,3)); B = np.zeros((3,3)); D = np.zeros((3,3))
    for k, angle in enumerate(angles):
        if angle is not None:
            for i in range(3):
                for j in range(3):
                    A[i][j] += Q_matrix[angle][i][j] * (z[k+1] - z[k])
                    B[i][j] += 1/2 * Q_matrix[angle][i][j] * (z[k+1]**2 - z[k]**2)
                    D[i][j] += 1/3 * Q_matrix[angle][i][j] * (z[k+1]**3 - z[k]**3)
    AB = np.hstack((A, B)); BD = np.hstack((B, D))
    ABD = np.vstack((AB, BD))
    return D, ABD, h

def getRandomLayup():
    return np.random.choice([0, 15, -15, 30,-30,45,-45,60,-60,75,-75,90,-90], size=(1, 4))

def getRandomPly():
    return np.random.choice([0, 15, -15, 30,-30,45,-45,60,-60,75,-75,90,-90])

def shuffleLayup():
    layup_ratio = [0,0, 90, 90, 15, -15, 45, -45, 45, -45, 45, -45, 45, -45, 60, -60, 60, -60, 75, -75, 30, -30, 30, -30, 30, -30]
    randomOrder = random.shuffle(layup_ratio)
    return layup_ratio

def returnMainLayup():
    return [45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 45, -45, 30, -30, 30, -30, 45, -45, 60, -60, 60, -60, 15, -15,  0, 90, 0, 90,45, -45, 90, 0]

def main():
    safe_ns = []
    base_layups = []
    bucklingStresses = []
    for i in range(0, 1, 1):
        #Works for getting ratio
        '''base_layup = getRandomLayup()[0]
        base_layup = base_layup.tolist()
        base_layups.append(base_layup)'''
        #layup = [getRandomPly()]
        """Start loop for getting the n at which it doesn't buckle"""
        #base_layup = shuffleLayup()
        #base_layup = [45, -45, 45, -45, 45, -45, 30, -30, 60, -60, 45, -45, 45, -45, 0, 90, 30, -30, 60, -60, 45, -45, 30, -30, 30, -30, 0, 90]
        base_layup = [45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 45, -45, 30, -30, 30, -30, 45, -45, 60, -60, 60, -60, 15, -15,  0, 90, 0, 90,45, -45, 90, 0]
        #base_layup = [45, -45, 90, 0]
        base_layups.append(base_layup)
        this_failed = True;
        n_ = 1
        while this_failed:
            """Make ABD Matrix"""
            D, ABD, h = getABDMatrix(base_layup, n_)
            
            d = np.linalg.inv(D)
            E1b = 12 / (h**3 * d[0][0])
            #print(f"E1b = {E1b}")
            #print(h)                                                        #mm
    
            a = 500
            b = 450 # can change these values around to see differences
            #k = 0
            #n = 2
            m = 1
            
            
            Iz = np.pi*(R**4 - (R-h)**4)/4
            print(f"Iz = {Iz}")
            actualStress = Mz * R / Iz  # N/(mm^2)
            
            print(f"Stress from moment = {actualStress}")

            
            '''If you want those old comments, insert them here'''
            
            
            bucklingStress = getBucklingLoadCCSS(D, m, a, b) / h  # N/mm / mm = N / (mm^2)
            
            print(f'Buckling Stress = {bucklingStress}')
            print(f'Buckling Load = {bucklingStress * h}')
            if bucklingStress > (1.5 * (actualStress)): #2x for factor of safety
            
                FOS = bucklingStress / actualStress
                this_failed = False;
                #len(layup)
                #safe_ns.append(len(layup))
                safe_ns.append(n_)
                bucklingStresses.append(bucklingStress)
            n_ += 1
            #layup.append(getRandomPly())
            #print(layup)
            

    min_n = np.min(safe_ns)
    min_weight_layups_w_load = [[base_layup, buckStress] for base_layup, n__, buckStress in zip(base_layups, safe_ns, bucklingStresses) if n__ == min_n]
    #print(min_weight_layups_w_load)
    #print(values_at_indexes)
    ordered_lightest_layups = sorted(min_weight_layups_w_load, key=lambda x: x[1])
    print(ordered_lightest_layups)
    print(f"Min N = {min_n}")
    print(f"num mins = {len(ordered_lightest_layups)}")
    numPlys = len(ordered_lightest_layups[-1][0])*2*min_n
    thicknessLaminate = numPlys * t
    weight = thicknessLaminate / 1000 * b/1000 * a/1000 * 1610
    print(f"Num plys = {numPlys}")
    print(f"FOS = {FOS}")
    print(f"Thickness = {thicknessLaminate} mm")
    print(f"Weight = {weight} kg")
    thicknessAluminium = 7 #mm
    weightOfAluminium = thicknessAluminium / 1000 * b/1000 * a/1000 * 2770
    print(f"Weight of Al = {weightOfAluminium} kg")
    
    #print(f"Num low weights = {np.size(ordered_lightest_layups)}")
    # Create a dictionary to store the counts of each value
    value_counts = {}
    #1, 1, 1, 3, 3
    #[0, 60, 15, 45, 45, 45, 30, 30, 30, -30, -30, -30, -45, -45, -45, -15, -60, 0]
    #take random ordering of above such that bucklingload is maximised. And that n is minimised where n is symetrising?

    #or does probabilility naturally do iit for you if you just add single plys until it's good?
    #also whats the bluckling case?
    
    #print(selected_values)
    
    #Value to beat = 0.16191011754224063
    
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










