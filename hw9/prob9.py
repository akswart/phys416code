# -*- coding: utf-8 -*-
"""
Created on Tue May  5 13:04:04 2020

@author: akswa
"""
import numpy as np

N = 5
c = 1


B = np.zeros((N,N))
C = np.zeros((N,N))


#Construct B diagonals
b1 = np.diag(np.ones(N-1),k=1)
b2 = np.diag(-np.ones(N-1),k=-1)
B += b1 +b2
# B corners
B[-1,0] = 1
B[0,-1] = -1

# Construct C
C = abs(B)

# IMPORTANT
# Note that for the Lax scheme, 
# the cfl condition is satisfied iff abs(c*tau/h) <= 1
# according to the textbook on pg 221


# Construct A
# pass it cfl = c*tau*h
def make_A(cfl):
    A = .5*C - (cfl/2)*B
    return A





###################################
# Part A
###################################
from power1 import power1

A1 = make_A(1.0001)
x=np.array([[i for i in range(1,N+1)]]).T

eigval, eigvec = power1(A1,x,1.0e-3,20)

#eigenvalue,eigenvector = np.linalg.eig(A1)
# get max eigenvalue
#emax = np.max(eigenvalue)
#emax_index = np.argmax(eigenvalue)
#evmax = eigenvector[:,emax_index]

#print(emax)
print(eigval)