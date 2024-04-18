#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 19:48:15 2024

Description:
# -- There are 42 total sectors of radians 0.15 in this fuselage. Each sector has 450mm width
This file splits the fuselage up into sectors and analyses the min viable thickness of plys
for that section. It reiterates over the fuselage several times, recomputing Iz and thickness if changed.
After running, the print shows the final weight of the fuselage for part B


@author: nathanrawiri
"""


"""Assignment Values"""
import Material_Constants as cons
E1, E2, G12, v12, t, Xt, Xc, Yt, Yc, S, v21, Q_basic = cons.main() #Get material values

import numpy as np
import Build_ABD_Matrix as abdFile
#import Final_Bottom_Case_PartB as bottom
import matplotlib.pyplot as plt
import Failure_Indexes as FI_file

rad_sec = 0.15                          # [rad] radians of a sector of 450mm width
R = 3000                                # (mm)  Radius of Fuselage

#Loads
M = 15e9                                # MNmm
V = 1.5e6                               # MN

def getABD(layup):
    D, ABD, _ = abdFile.getABDMatrix(layup, 1)
    return D, ABD

def getOriginalBottomLayup(): # for now can use this original layup. May change later to include bottom section in loop.
    #return [45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 90, 0, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30\
        #, -30, 45, -45, 30, -30, 30, -30, 45, -45, 60, -60, 60, -60, 0, 90, 0, 15, -15, 0, 90, 0, 90, 45, -45, 90, 0]
        # The commented one above was the assumed bottom layup before integrating this strong loop. Now we get the below because Iz is more accurate defined in this loop.
    #return [45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 30, -30, 30, -30, 30, -30, 60, -60, 60, -60, 60, -60, 60, -60, 15, -15, 0, 90, 75, -75, 0, 90, 0, 90]
                 #, 45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 30, -30, 30, -30, 60, -60, 60, -60, 60, -60, 15, -15, 15, -15, 0, 90, 75, -75, 0, 90, 0, 90, 0, 90]
     return [45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 30, -30, 30, -30, 30, -30, 60, -60, 60, -60, 60, -60, 60, -60, 15, -15, 0, 90, 75, -75, 0, 90, 0, 90\
                  , 30, -30, 30, -30, 60, -60, 60, -60,  15, -15, 0, 90, 0, 90, 45, -45, 0, 90]

    
def getTopLayup(): #same as above. May change later
    return [0, 0, -45, 45, 0, 0, 90]

def Iz_sec(t, alpha):
    theta = np.pi/4 - (alpha - rad_sec); phi = np.pi/4 - alpha;
    return ((R**4 - (R-t)**4)/4) * ((theta - 0.5*np.sin(2*theta))-(phi - 0.5 * np.sin(2 * phi))) /2

def y_(t, alph):
    return (2*np.sin(alph) / (3 * alph)) * ((R**3-(R-t)**3)/(R**2-(R-t)**2))

def A_(alph, t):
    return alph * (R**2 - (R-t)**2)
    
def Tau_(I, t, alph):
    b = 2*t; 
    Q = y_(t, alph) * A_(alph, t);
    return V * Q / (I*b)

def M_(I, y):
    return M * y / I

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

def getBucklingLoadShear(D, a, b):
    D11 = D[0][0]; D12 = D[0][1]; D66 = D[2][2]; D22 = D[1][1];
    return 9 * np.pi**4 * b / (32 * a**3 * (D11 + 2*D12 + 2*D66)*a**2/b**2 + D22 * a**4 / b**4) / 1.27

def getBucklingLoadCompAndShear(D, a, b, k):
    D11 = D[0][0]; D12 = D[0][1]; D66 = D[2][2]; D22 = D[1][1];
    plus = (np.pi**2/a**2) * ((D11 + 2 * (D12+2*D66)*(a**2/b**2) + D22*a**4/b**4)\
                /(2 - 8192/81*a**2*k**2/(b**2*np.pi**4))) * 5 + np.sqrt(9+(65536/81)*(a**2/(b**2*np.pi**4)*k**2))
    minus = (np.pi**2/a**2) * ((D11 + 2 * (D12+2*D66)*(a**2/b**2) + D22*a**4/b**4)\
                /(2 - 8192/81*a**2*k**2/(b**2*np.pi**4))) * 5 - np.sqrt(9+(65536/81)*(a**2/(b**2*np.pi**4)*k**2))
    return np.max((plus, minus))

"""Initialise Iz, t, layup arrays"""
# 21 is half of the way around the circle. But adding an extra because top and bottom are wider circles
# potential point of error if 20 is not enough sections.
# But better to do it this way than appending because this way is more scalable for re-looping
Iz = [2.1e12] * 22 #starting with what Iz / 42 if Iz is for constant t = 2.5mm all around
ts = [2.5] * 22 # Starting with 2.5mm average thickness all around
ts[0] = 9.72
layups = [[] for _ in range(22)] 
len_layups = [[] for _ in range(22)] 
weight_section = [[] for _ in range(22)] 


layups[0] = getOriginalBottomLayup()    # make bottom layup 1st layup
layups[21] = getTopLayup()              # make top layup last layup
Iz_total = Iz[0] + Iz[-1] + sum(Iz[1:-1])*2  # Total Iz of the fuselage. Subject to change in loop
len_layups[0] = len(layups[0]); len_layups[21] = len(layups[21]); 


def increaseLayup(i):
    layups[i] = [45, -45] + layups[i]
    #print(f"Length = {len(layups[i])}")
    ts[i] = len(layups[i]) * 0.135 * 2   # 2 because symmetrical
    Iz[i] = Iz_sec(ts[i], alphas[i])
    
def decreaseLayup(i):
    layups[i] = layups[i][2:]           # this removes the 1st 2 elements from the list
    

Taus = []; M_s = []; ys = [];           # Declare arrays for stress distribution plots
Nss = []; Nxs = []                      # Declare arrays for force distribution plots

alphas = []                             # Make list of the max angle for each segement
alpha = 0.15/2 + 0.15                   #the first alpha
while alpha < np.pi:
    alphas.append(alpha); alpha += 0.15

Iz_totals = []
#while Iz_totals[-1] - Iz_totals[-2] < 5e11: # change this to a while loop checking tolerances
for i in range(1, 7, 1):
    i = 1
    Taus = []; M_s = []; ys = []; Nss = []; Nxs = []; Ns_ratios = []; combined_bucklings = []; comps = []
    for alpha in alphas:
        layups[i] = layups[i-1]                 # start by making this one the same as the previous
        y = -R*np.cos(alpha)                    # get the y value at the sector.
        """Start loop to adjust layup of sector"""
        adjusting_layup = True;                 # Set True until layup for sector confirmed
        while adjusting_layup:                  # Until layup confirmed, repeat loop
            #print(y)
            #print(alpha)                       # uncomment these prints to ensure in the right ranges
            Tau_sec = Tau_(Iz_total*4, ts[i], alpha) # Get shear stress at section
            Mz = M_(Iz_total*4, y)              # Get Stress due to moment at section
            Ns = Tau_sec * ts[i]; Nx = Mz * ts[i]       # Get Axial loads in  N/mm.
            D, ABD = getABD(layups[i])          # Get ABD Matrix
            """Check Failure Indexes"""
            FI = FI_file.checkFI(Nx, Ns, ABD, layups[i]) #Get CLT Failure Index utilising Puck and Max. Stress Failure Criteria
            m = 1; a = 500; b = 450; k  = 1                    # Setting width and length of arbitrary plate approximation
            bucklingLoad = getBucklingLoadCCSS(D, 1, 500, 450)  # Buckling load in N / mm 
            shearBL = getBucklingLoadShear(D, 500, 450)         # Shear buckling load in N/ mm
            shearFOS = shearBL / abs(Ns)
            Nx_ratio = Nx / bucklingLoad
            Ns_ratio = (1/abs(Nx_ratio))**0.5
            combined_buckling = getBucklingLoadCompAndShear(D, a, b, k)
            FI_FOS = 1/FI
            """Several safety checks below. Either increases or decreases layup based on SFs"""
            if FI_FOS < 2:
                increaseLayup(i)
            elif y < 0 and abs(combined_buckling) < 2 * (abs(Nx)+Ns):  # In bottom half, check combined buckling load
                increaseLayup(i)
            elif y < 0 and bucklingLoad > 2.3 * abs(Nx) and len(layups[i])>4:
                decreaseLayup(i)
            elif y < -1 and bucklingLoad < 1.5 * abs(Nx):
                increaseLayup(i)
            elif y > 0 and (abs(combined_buckling) < 3 * (Nx+Ns)): # and FI < 0.33: #and FI > 0.33:
                increaseLayup(i)
                print(f"buck: {abs(combined_buckling)}: Nx: {Nx}, Ns: {Ns}")
            elif y > 0 and abs(combined_buckling) > 3 * (Nx+Ns) and len(layups[i])>4 and FI < 0.33: #and FI > 0.33:
                decreaseLayup(i)
            elif y > 0 and FI_FOS < 400:
                increaseLayup(i)
            else:
                adjusting_layup = False         # FI's are all good. Layup is good. Move to next sector
        len_layups[i] = len(layups[i])
        ys.append(y)                            # append for stress distribution curve 
        Taus.append(Tau_sec);M_s.append(Mz);    # append for stress distribution curve 
        Nxs.append(Nx); Nss.append(Ns)          # append for force distribution curve
        #print(f"Shear FOS = {shearFOS}")
        #print(f"Buckling of Combined Compresssion and shear = {combined_buckling}")
        #print(f"Buckling of Compresssion = {bucklingLoad}")
        #print(f"Actual load = {Nx / t * ts[i]}")
        combined_bucklings.append(combined_buckling)
        comps.append(bucklingLoad)
        Ns_ratios.append(Ns_ratio)
        ts[i] = len(layups[i]) * t * 2          # thickness of the laminate * 2 because symmetrical
        weight_section[i] = len_layups[i] * 2 * 0.135 * 0.5 * 0.45 / 1000 * 1610 
        #Finished this iteration of loop
        alpha += 0.15; i += 1                   # move to the next sector
    """After looping thorugh the fuselage"""
    Iz[0] = Iz[1]; Iz[-1] = Iz[-2]
    Iz_total = Iz[0] + Iz[-1] + sum(Iz[1:-1])*2 # *2 because both sides of fuse
    ts[0] = ts[1]; ts[-1] = ts[-2]
    Iz_totals.append(Iz_total)
#print(f"Iz_final = {Iz_total}")
    
weight_section[0] = weight_section[1]; weight_section[-1] = weight_section[-2]
weight_total = weight_section[0] + weight_section[-1] + sum(weight_section[1:-1])*2
print(f"Weight Composite = {weight_total:.2f} kg")
weight_Al = 2 * np.pi * 3 * 7/1000 * 0.5 * 2770 
print(f"Weight aluminium = {weight_Al:.2f} kg")
print(f'Weight reduction = {100-(100*weight_total/weight_Al):.2f}%')

print(f"Thickness per section from bottom to top in radians 0.15 increments: {len_layups}")



plt.figure()
plt.plot(Nss, np.array(ys)/1000, label = "Ns", color='orange')
plt.plot(Nxs, np.array(ys)/1000, label = "Nx", color='blue')
plt.xlabel('Applied Force (N/mm)')
plt.ylabel('y (m)')
#plt.title('Freq of angle present in lowest FI layups for tensile Max')
plt.grid(True)
#plt.ylim(ymin=0)#, ymax=150)  # Adjust the limits for the y axis
plt.legend()
plt.show()

plt.figure()
plt.plot(Taus, np.array(ys)/1000, label = "Tau", color='orange')
plt.plot(M_s, np.array(ys)/1000, label = "Moment_Stress", color='blue')
plt.xlabel('Stress (MPa)')
plt.ylabel('y (m)')
#plt.title('Freq of angle present in lowest FI layups for tensile Max')
plt.grid(True)
#plt.ylim(ymin=0)#, ymax=150)  # Adjust the limits for the y axis
plt.legend()
plt.show()





