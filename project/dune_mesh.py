#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 22:10:19 2020

@author: aswart

Generates mesh for navier_stokes_dune.py
For whever reason, this is really sensitive to Seg faults when tweaking the 
parameters on the tails of the gaussian. Not sure why.

"""

from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt



def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def gaussian_mesh(x0,x1,y1,mu,sig,mesh_size=24):
    x_val = np.linspace(x0,x1,mesh_size) #100 points is probably more than enough
    g = y1*gaussian(x_val,mu,sig)
    g[np.abs(g) < 10**-4] = 0  # Set small values to zero
    
    gaussian_mesh_vertices = [Point(x0,0.0)] # Append beg point
    for i in range(1,len(g)-1):
        gaussian_mesh_vertices.append(Point(x_val[i],g[i]))
    gaussian_mesh_vertices.append(Point(x_val[-1],0.0)) # Drop down to zero
    gaussian_mesh_vertices.append(Point(x0,0.0)) # And loop back to beginning for closed manifold
    gaussian_mesh_vertices.reverse() #Points need to be done in counter-clockwise order
    
    return gaussian_mesh_vertices

def generate_dune_mesh(large_scale = True,mesh_size=32):
    # Keeping dune size the same, we zoom in or out
    if large_scale: # Zoom out x=(-4,20) y=(0,10)
        channel = Rectangle(Point(0, 0), Point(2.2, .41))
    else: # Zoom in x=(-4,10) y = (0,2)
        channel = Rectangle(Point(-4, 0), Point(10,2))
    dune_mesh = gaussian_mesh(.1,.7,.12,.45,.1,int(.75*mesh_size))
    dune = Polygon(dune_mesh)
    gmesh = generate_mesh(channel-dune,mesh_size)
    plot(gmesh)
    return gmesh

if __name__ == "__main__":
    dune_mesh = generate_dune_mesh(True,32)
    plot(dune_mesh)






