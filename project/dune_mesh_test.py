#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 22:10:19 2020

@author: aswart

just messing around with meshes
"""

from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt



def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def gaussian_mesh(x0,x1,y1,mu,sig):
    x_val = np.linspace(x0,x1,50) #100 points is probably more than enough
    g = y1*gaussian(x_val,mu,sig)
    g[np.abs(g) < 10**-4] = 0  # Set small values to zero
    gaussian_mesh_vertices = [Point(x0,0.0)] # Append beg point
    for i in range(len(g)):
        gaussian_mesh_vertices.append(Point(x_val[i],g[i]))
    gaussian_mesh_vertices.append(Point(x_val[-1],0.0)) # Drop down to zero
    gaussian_mesh_vertices.append(Point(x0,0.0)) # And loop back to beginning for closed manifold
    gaussian_mesh_vertices.reverse()
    return gaussian_mesh_vertices



'''
domain_vertices = [Point(0.0, 0.0),
                   Point(0.0, 0.0),
                   Point(10.0, 0.0),
                   Point(10.0, 2.0),
                   Point(8.0, 2.0),
                   Point(7.5, 1.0),
                   Point(2.5, 1.0),
                   Point(2.0, 4.0),
                   Point(0.0, 4.0),
                   Point(0.0, 0.0)]


channel = Rectangle(Point(-4, 0), Point(20, 10))
#dune = Polygon(domain_vertices)
#mesh = generate_mesh(channel-dune,64)
#plot(mesh)
'''

channel = Rectangle(Point(-4, 0), Point(20, 10))

dune_mesh = gaussian_mesh(0,10,1,2,1)
dune = Polygon(dune_mesh)
gmesh = generate_mesh(channel-dune,32)
plot(gmesh)



