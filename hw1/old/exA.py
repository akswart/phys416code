# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:40:52 2020

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
    while num_iter < cutoff:
        num_iter += 1
        n = hailstone(n)
        hailstone_list.append(n)
        print(n)
        
        if n == 1:
            return hailstone_list
    return hailstone_list
    


if __name__ == '__main__':
    h = hailstone_iterator(27)
    x = range(len(h))
    
    plt.plot(x,h)
    print(h)