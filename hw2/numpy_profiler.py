# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 14:58:43 2020

@author: akswa
"""
import numpy as np


# Test speed of numpy append
def np_append():
#    b = np.array([])
    for j in range(1000):
        a = np.array([-1])
        for i in range(1000):
            a = np.append(a,i)
        if j == 0:
            b = a
        else:
            b = np.vstack((b,a))
    return b        

# Test speed of python append then casting to numpy array
# Test speed of numpy append
def py_append():
    b = []
    for j in range(1000):
        a = [-1]
        for i in range(1000):
            a.append(i)
        b.append(a)
    return np.array(b) 


if __name__ == '__main__':
    import timeit
    '''
    a = np_append()
    b = py_append()
    print(a.shape,b.shape)
    print(np.all(a==b))
    
    '''
    #"""
    t1 = timeit.timeit("np_append()",setup="from numpy_profiler import np_append",number = 10)
    t2 = timeit.timeit("py_append()",setup="from numpy_profiler import py_append",number = 10)
    print('Np',t1)
    print('Py',t2)
    print(t1/t2)
    
    #"""