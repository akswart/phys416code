# python 3 version 2/15
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
def force_sphere(h,rho):
# force term for the floating spshere problem'
# frt 2/14
    f = 4*rho - (h**2)*(3-h);
    return f

rho=float(input(' please input rho ( < 1): '))
assert rho<=1 and rho>=0.0 ,'Rho must be between 0 and 1'

h_start = 1
hold_start = 0 # for the secant method

h = h_start; # first guess
eps = 0.0
error_sol_newton =  np.abs(force_sphere(h,rho))
error_iter_new = np.array([error_sol_newton])
error = 1e-7 # set error tolerance

itermax = 1000 # max number of iterations

for iter in range(1,itermax):
    if(eps > error or error_sol_newton > error):
        hnew = h - force_sphere(h,rho)/(3*h**2-6*h)
        eps = np.abs(hnew - h)
        error_sol_newton = np.abs(force_sphere(hnew,rho))
        h = hnew
        error_iter_new = np.append(error_iter_new,error_sol_newton)
    else:
        break

hnewton=h

print('Newton, solution at iteration number %f the error is %f and the depth = %f '%(iter,eps,h))

# -----------now solve the same equation using the secant method---------

h = h_start; # first guess
h_old = hold_start; # need an hold

eps = 0;
error_sol_secant = np.abs(force_sphere(h,rho));
error_iter_sec = np.array([error_sol_secant])

itermax = 1000; # max number of iterations

for iter in range(1,itermax):
    if(eps > error or error_sol_secant > error):
        hnew = h - force_sphere(h,rho)/(force_sphere(h,rho)-force_sphere(h_old,rho))*(h-h_old);
        eps = abs(hnew - h)
        error_sol_secant = abs(force_sphere(hnew,rho))
        h = hnew;
        error_iter_sec = np.append(error_iter_sec,error_sol_secant)
    else:
        break


print('Secant, solution at iteration number %f the error is %f and the depth = %f '%(iter,eps,h))

# now plot the solition

hplot = np.arange(0,2*h,0.1*h)

plt.figure(1)
plt.plot(hplot,force_sphere(hplot,rho),'-')
plt.plot(hnewton,force_sphere(hnewton,rho),'+')
plt.xlabel('depth')
plt.ylabel('force curve')
plt.legend(['force curve','solution'])
plt.title('solution from Newton')
plt.grid(True)

# generate a circle, radius 1 centered at  (1-h)
plt.figure(2)
# anglerange=np.arange(0,2*np.pi,np.pi/20)
anglerange=np.linspace(0,2*np.pi,endpoint=True)
plt.plot(np.sin(anglerange),1-hnewton+np.cos(anglerange),'g-')
plt.plot(np.arange(-2,0,2),'b-')
plt.title('solution from Newton')
plt.axis('equal')
plt.grid(True)


itermax=np.max([len(error_iter_new),len(error_iter_sec)])
print(itermax)
# % now plot
if(itermax > 1):
    plt.figure(3)
    plt.semilogy(np.arange(0,len(error_iter_new)),error_iter_new,'r-')
    plt.semilogy(np.arange(0,len(error_iter_sec)),error_iter_sec,'b-')
    plt.xlabel('iteration')
    plt.ylabel('error')
    plt.legend(['Newton','Secant'])
    plt.axis([1,itermax,0, 2])
    plt.grid(True)
else:
    print('Lucky guess, you got it first time!')


plt.show()
