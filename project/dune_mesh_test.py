#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 22:10:19 2020

@author: aswart

just messing around with meshes
"""

from fenics import *
from mshr import *
mesh = Mesh()
domain_vertices = [Point(0.0, 0.0),
                   Point(10.0, 0.0),
                   Point(10.0, 2.0),
                   Point(8.0, 2.0),
                   Point(7.5, 1.0),
                   Point(2.5, 1.0),
                   Point(2.0, 4.0),
                   Point(0.0, 4.0),
                   Point(0.0, 0.0)]

channel = Rectangle(Point(-4, 0), Point(20, 10))
dune = Polygon(domain_vertices)
mesh = generate_mesh(channel-dune,40)
plot(mesh)
