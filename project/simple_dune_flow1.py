# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 16:58:13 2020

@author: akswa

Objective is to create a simple fluid flow simulation over a gaussian curve

"""


from fluidsim.solvers.ns2d.solver import Simul
params = Simul.create_default_params()
[attr for attr in dir(params) if not attr.startswith('_')]




params.nu_2 = 1e-3
params.forcing.enable = False

params.init_fields.type = 'noise'

params.output.periods_save.spatial_means = 1.
params.output.periods_save.spectra = 1.
params.output.periods_save.phys_fields = 2.

sim = Simul(params)