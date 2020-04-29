#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 22:06:54 2020

@author: aswart

Graphing the output of hw8_exC_parallel_step_opt.py given in outputfile.csv

The provided output.csv is from a run done on my server with the following parameters:
    beg_scale = .1
    end_scale = 5
    num_point = 500

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.close('all')
df = pd.read_csv("outputfile.csv")
#plt.figure()
df.plot(kind='scatter',x='Scale',y='Time')
plt.xlabel("Multiplier from suggested values for grid size")
#plt.figure()
df.plot(kind='scatter',x='Scale',y='Iterations')
plt.xlabel("Multiplier from suggested values for grid size")