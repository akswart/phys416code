# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:13:29 2020

@author: akswa
"""

import numpy as np
import matplotlib.pyplot as plt

def hailstone(n):
    if n%2 == 0:
        return int(n/2)
    elif n%2 == 1:
        return 3*n+1
    elif type(n) != int:
        raise TypeError("Input is not integer")

def hailstone_iterator(n):
    num_iter = 0
    hailstone_list = [n]
    cutoff = 10000
    while num_iter < cutoff and n != 1:
        num_iter += 1
        n = hailstone(n)
        hailstone_list.append(n)

    return hailstone_list
    
def vectorized_hailstone(N):
    s = []
    f = []
    g = []
    for n in range(1,N+1):
        h = hailstone_iterator(n)
        s += (np.ones(len(h),dtype=int)*n).tolist()
        f += h
        g.append(len(h)-1)
    return s,f,g

if __name__ == '__main__':
    N = 200
    s,f,g = vectorized_hailstone(N)
    fig, ax = plt.subplots(1,2,gridspec_kw={'width_ratios': [2, 1]})
    ax[0].scatter(s,f,s=1,c='C0')
    ax[1].scatter(range(1,N+1),g,s=2,color='red')
    
    ax[0].set_title('range of computed values during iteration')
    ax[1].set_title('number of iterations')
    
    for i in ax:
        i.set_xlim(1,N)
        i.set_ylim(1,N)
        