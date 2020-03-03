# square root finding algorithm
import numpy as np
import matplotlib.pyplot as plt

xr=float(input('input number: '))

x0 = xr -1
x = x0 # first guess

error_min = 1.0e-6
error = 100

counter = 0
counter_max = 100
xplot=np.array([])

while( error > error_min and counter < counter_max):
    counter = counter + 1;

    xplot = np.append(xplot,x)

    dx = (x**2-xr)/(2*x)
    x = x - dx

    error = abs(x**2-xr);


if (counter == counter_max):
    print(' Max counts exceeded, no solution found')
else:
    print('The square root of ',xr,' is: ',x);

plt.figure(1)
plt.plot(xplot)
plt.xlabel('iterations')
plt.ylabel('solution')
plt.show()