# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt

# 4 Masses, 4 forces
#F1,F2,F3,F4) = (0,0,0,0)

# Initial pos
#(x1,x2,x3,x4) = (1,2,3,4)

# Sign convention assumes springs are in tension
# i.e. the springs pull in, wanting to be shorter
# Forces to left are neg, forces to right are positive

# Set spring constants
ratios = np.logspace(-3,3,100)
#k1 = 1
lengths = []
for k1 in ratios:
    k2 = 1
    # Rest length is 1 for all springs
    L = 1
    
    #F1 = k1*(x2-x1-L) + k2*(x3-x1-L)
    #F2 = -k1*(x2-x1-L) + k1*(x3-x2-L) + k2*(x4-x2-L)
    #F3 = -k1*(x3-x2-L) - k2*(x3-x1-L) + k1*(x4-x3-L)
    #F4 = -k1*(x4-x3-L) - k2*(x4-x2-L)
    
    # Set arrays
    #K_old = np.array([ [-k1-k2 ,k1       ,k2       ,0],
    #               [k1     ,-2*k1-k2 ,k1       ,k2],
    #               [k2     ,k1       ,-2*k1-k2 ,k1],
    #               [0      ,k2       ,k1       ,-k1-k2] ])
    
    #b_old = np.array([[-k1*L-k2*L,-k2*L,k2*L,k1*L+k2*L]]).T
    # Since K_old is singular, we need to fix a varibale, fix x1 at 0
    
    
    K = np.array([ [1 ,0       ,0       ,0],
                   [k1     ,-2*k1-k2 ,k1       ,k2],
                   [k2     ,k1       ,-2*k1-k2 ,k1],
                   [0      ,k2       ,k1       ,-k1-k2] ])
    
    b = np.array([[0,-k2*L,k2*L,k1*L+k2*L]]).T
    
    #print(np.linalg.inv(K))
    # Solve Kx = b
    x1 = np.linalg.inv(K) @ -b #K^-1 * b
    #x2 = np.linalg.solve(K,b)    
    lengths.append(x1[-1])
    #print(x2)

plt.semilogx(ratios,lengths)

