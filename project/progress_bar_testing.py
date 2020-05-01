#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 12:44:13 2020

@author: aswart
 messing around with progress bars


"""

import time
from tqdm import tqdm

for i in tqdm(range(20)):
    time.sleep(.5)
    
