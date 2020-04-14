# neutrn2 - Slightly modified version of neutrn.m program 
# to solve the neutron diffusion equation 
# using the Forward Time Centered Space (FTCS) scheme.
# clear help neutrn2  # Clear memory and print header
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

#* Initialize parameters (time step, grid points, etc.).

N = int(input('Enter the number of grid points: '))

D = 1.   # Diffusion coefficient
C = 1.   # Generation rate
L_critical = np.pi*np.sqrt(D/C)
print(' Critical length is ',L_critical)

L = float(input('Enter system length: '))
# The system extends from x=-L/2 to x=L/2
h = L/(N-1.)  # Grid size

tau_max = 0.5*h**2/D
print(' Maximum timestep = ',tau_max)
tau = float(input('Enter time step: '))

coeff = D*tau/h**2
coeff2 = C*tau
if( coeff < 0.5 ):
  print('Solution is expected to be stable')
else:
  print('WARNING: Solution is expected to be unstable')

#* Set initial and boundary conditions.
nn = np.zeros(N)        # Initialize density to zero at all points
nn_new = np.zeros(N)    # Initialize temporary array used by FTCS
nn[round((N-1)/2)] = 1./h   # Initial cond. is delta function in center
## The boundary conditions are nn(1) = nn(N) = 0

#* Set up loop and plot variables.
xplot = np.arange(0,N)*h-L/2.   # Record the x scale for plots
iplot = 1                 # Counter used to count plots
nstep = int(input('Enter number of time steps: '))
nplots = 50               # Number of snapshots (plots) to take
plot_step = nstep/nplots  # Number of time steps between plots

nnplot = np.array([])
tplot = np.array([])
nAve = np.array([])

#* Loop over the desired number of time steps.
for istep in range(0,nstep+1):  ## MAIN LOOP ##

  #* Compute the new density using FTCS scheme.
  nn_new[1:(N-1)] = nn[1:(N-1)] + coeff*(nn[2:N] + nn[0:(N-2)] - 2*nn[1:(N-1)]) + coeff2*nn[1:(N-1)]
  nn = nn_new        # Reset temperature to new values
  
  #* Periodically record the density for plotting.
  if( istep%plot_step < 1 or istep ==0 ):   # Every plot_step steps
    if(iplot==1):
        nnplot = np.copy(nn)
    else:
        nnplot = np.vstack((nnplot,nn))       # record nn(i) for plotting
    tplot =  np.append(tplot,istep*tau)     # Record time for plots
    nAve= np.append(nAve,np.mean(nn))        # Record average density
    iplot = iplot+1
    print('Finished %d of %d steps'%(istep,nstep))

tt,xx = np.meshgrid(xplot,tplot)
#* Plot density versus x and t as a 3D-surface plot
fig = plt.figure(1)
plt.clf()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xx, tt, nnplot, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
ax.set_xlabel('Time')
ax.set_ylabel('x')
ax.set_zlabel('n(x,t)')
# title('Neutron diffusion')

# compute alpha1, which is the growth rate of the fastest mode
# and overplot with the average coming out of the code
alpha1 = C -D *np.pi**2/L**2
nAve_theory = nAve[0]*np.exp(alpha1*tplot)
#* Plot average neutron density versus time
plt.figure(2)
plt.clf()
plt.plot(tplot,nAve,'*',tplot,nAve_theory,'-')
plt.xlabel('Time')
plt.ylabel('Average density')
plt.title('L = '+str(L)+r'( $L_c = \pi$ )')
plt.legend(['model',r'theory using $\alpha _1$'])
plt.show()
