#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 22:22:16 2024

@author: nathanrawiri
"""

E1 = 140
E2 = 11.2  
G12 = 5
v12 = 0.33

t = 0.135
#Need strength
Xt = 2200 # MPa
Xc = 1800 #MPa
Yt = 70 #MPa
Yc = 300 #MPa
S = 100 #MPa

v21 = v12 * E2/E1
Q = 1-v12*v21
Q11 = E1 / Q 
Q22 = E2 / Q
Q12 = v12*E2 / Q
Q66 = G12 
Q_basic = [[Q11, Q12, 0], [Q12, Q22, 0],[0, 0, Q66]]


def main():
    return E1, E2, G12, v12, t, Xt, Xc, Yt, Yc, S, v21, Q_basic