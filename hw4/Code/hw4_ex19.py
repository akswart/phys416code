#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:15:57 2020

@author: aswart
"""
import numpy as np

#SEE SCANNED PAPER FOR CALCULATIONS AND MORE EXPLANATION

# Use newtons method to find maximum
# Want to solve for lambda where du/dL = 0
# Let all constants equal 1
c = 1
h = 1
T = 1
k = 1

g = lambda x: x - np.exp(1/x)/(5*(np.exp(1/x)-1))
g_p = lambda x: 1 - 1/(5*x**2*(np.exp(1/x)-1)) - 1/(5*x**2*(np.exp(1/x)-1)**2)

# Time for newtons method
tol = 10**-7
x0 = 1
x1 = .2
while abs(x0-x1) >= tol:
    x0 = x1
    x1 = x0 - g(x0)/g_p(x0)
    
print("Alpha is:",x1)
