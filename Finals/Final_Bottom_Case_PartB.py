#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 09:56 2024

Description: This file contains analyses of the chosen layup for the bottom  of the fuselage

@author: nathanrawiri


"""



import numpy as np
import Build_ABD_Matrix as abdFile

import Material_Constants as cons
E1, E2, G12, v12, t, Xt, Xc, Yt, Yc, S, v21, Q_basic = cons.main() #Get material values


"""Initialise loads"""

Ny = 0; Ns = 0; Mx = 0; My = 0; Ms = 0;
Mz = 15 * 10**9


R = 3000

def getBucklingLoadCCSS(D, m, a, b): #b is width, a is length
    D11 = D[0][0]; D12 = D[0][1]; D66 = D[2][2]; D22 = D[1][1]; lamb = (a/b)*(D22/D11)**(1/4)
    if (0 < lamb < 1.662):
        K = (4/lamb**2)+((2*(D12+2*D66))/(np.sqrt(D11*D22)))+((3/4)*lamb**2)
    elif (lamb >= 1.662):
        K = ((m**4+8*m**2+1)/(lamb**2*(m**2+1)))+((2*(D12+2*D66))/(np.sqrt(D11*D22)))+((lamb**2)/(m**2+1))
    return (np.pi**2/(b**2))*np.sqrt(D11*D22)*K

def getBucklingLoadSSSS(D, m, a, b):
    D11 = D[0][0]; D12 = D[0][1]; D66 = D[2][2]; D22 = D[1][1]; AR = a/b;
    return np.pi**2 * (D11*m**4 + 2*(D12+2*D66)*m**2*AR**2 + D22*AR**4 ) / (a**2 * m**2)
    
def getBucklingLoadCompAndShear(D, a, b, k):
    D11 = D[0][0]; D12 = D[0][1]; D66 = D[2][2]; D22 = D[1][1];
    plus = (np.pi**2/a**2) * ((D11 + 2 * (D12+2*D66)*(a**2/b**2) + D22*a**4/b**4)\
                /(2 - 8192/81*a**2*k**2/(b**2*np.pi**4))) * 5 + np.sqrt(9+(65536/81)*(a**2/(b**2*np.pi**4)*k**2))
    minus = (np.pi**2/a**2) * ((D11 + 2 * (D12+2*D66)*(a**2/b**2) + D22*a**4/b**4)\
                /(2 - 8192/81*a**2*k**2/(b**2*np.pi**4))) * 5 - np.sqrt(9+(65536/81)*(a**2/(b**2*np.pi**4)*k**2))
    return np.max((plus, minus))

"""MAIN FUNCTION"""
def main():

    #base_layup = [45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45,90, 0, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 45, -45, 30, -30, 30, -30, 45, -45, 60, -60, 60, -60, 0, 90, 0, 15, -15,  0, 90, 0, 90,45, -45, 90, 0]
    #base_layup= [45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 30, -30, 30, -30, 30, -30, 60, -60, 60, -60, 60, -60, 60, -60, 15, -15, 15, -15, 75, -75, 0, 90, 90, 0\
    #             ,45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 30, -30, 30, -30, 30, -30, 60, -60, 60, -60, 60, -60, 60, -60, 15, -15, 15, -15, 75, -75, 0, 90, 90, 0]
                # ,45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 30, -30, 30, -30, 30, -30, 60, -60, 60, -60, 60, -60, 60, -60, 15, -15, 15, -15, 75, -75, 0, 90, 90, 0]
                 #, 45, -45, 45, -45, 45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 30, -30, 45, -45, 45, -45, 45, -45, 60, -60, 60, -60, 60, -60, 15, -15, 75, -75, 0, 90, 90, 0, 0, 90, 90, 0]
    #base_layup = [45, -45, 90, 0]
    base_layup = [45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 30, -30, 30, -30, 30, -30, 60, -60, 60, -60, 60, -60, 60, -60, 15, -15, 0, 90, 75, -75, 0, 90, 0, 90\
                 , 45, -45, 45, -45,  30, -30, 30, -30, 60, -60, 60, -60,  15, -15, 0, 90, 0, 90, 0, 90]
    this_failed = True;
    n_ = 1
    while this_failed:
        """Make ABD Matrix"""
        D, ABD, h = abdFile.getABDMatrix(base_layup, n_)

        a = 500; b = 450; m = 1
        
        Iz = 210000000000000
        actualStress = Mz * R / Iz  # N/(mm^2)
        
        bucklingStress = getBucklingLoadCCSS(D, m, a, b) / h  # N/mm / mm = N / (mm^2)
        #both = getBucklingLoadCompAndShear(D, a, b, k)
        
        
        if bucklingStress > (2 * (actualStress)): #2x for factor of safety
            # This is safe.
            FOS = bucklingStress / actualStress
            this_failed = False;
            break
        else:
            FOS = bucklingStress / actualStress
            print(f"FOS = {FOS}")
            n_ += 1

    print(base_layup)
    print(f"Min N = {n_}")
    numPlys = len(base_layup)*2*n_
    thicknessLaminate = numPlys * t
    weight = thicknessLaminate / 1000 * b/1000 * a/1000 * 1610
    print(f"Num plys = {numPlys}")
    print(f"FOS = {FOS}")
    print(f"Thickness = {thicknessLaminate} mm")
    print(f"Weight = {weight} kg")
    thicknessAluminium = 7 #mm
    weightOfAluminium = thicknessAluminium / 1000 * b/1000 * a/1000 * 2770
    print(f"Weight of Al = {weightOfAluminium} kg")
    
    return 

main() 










