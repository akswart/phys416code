# -*- coding: utf-8 -*-
"""
Phys 416
HW 1
Alexander Swart
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

h = 10**-(np.arange(15,-.5,-.5))
delta = lambda h: np.abs(3 - ( (1+h)**3 - (1-h)**3 ) / (2*h) )

h_exact = 10**-(np.arange(15,-.05,-.05))
approx_error = lambda h: h**2 + sys.float_info.epsilon/h


fig = plt.figure()

plt.loglog(h,delta(h),marker='.',linestyle='',markersize=10)
plt.loglog(h_exact,approx_error(h_exact))

plt.title(r'error for the derivative of $x^{3}$ at x = 1.0')

plt.xlabel("h")
plt.ylabel(r'$\Delta(h)$')
plt.xticks([10**-14,10**-11,10**-8,10**-5,10**-2])
plt.yticks([10**-10,10**-8,10**-6,10**-4,10**-2,10**0])
plt.grid(True)
