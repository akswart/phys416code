from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
# advect - Program to solve the advection equation 
# using the various hyperbolic PDE schemes
# clear all  help advect  # Clear memory and print header
plt.close('all')
#* Select numerical parameters (time step, grid spacing, etc.).
method = 1
#N = int(input('Enter number of grid points: '))
N = 100
L = 1.     # System size
h = L/N    # Grid spacing
c = 1.      # Wave speed
print('Time for wave to move one grid spacing is ',h/c)
#tau = float(input('Enter time step: '))
tau = .2*h/c

coeff = (c*tau/h)**2    # Coefficient used by L-W scheme
print(coeff)
print('Wave circles system in %6.2f steps'%(L/(c*tau)))
#nStep = int(input('Enter number of steps: '))
nStep = int(.2*L/(c*tau))
#* Set initial and boundary conditions.
sigma = 0.1              # Width of the Gaussian pulse
k_wave = np.pi/sigma        # Wave number of the cosine
x = (np.arange(-1,N+1)+1./2.)*h - L/2  # Coordinates of grid points
# Initial condition is a Gaussian-cosine pulse
x0 = 0
a =  np.cos(k_wave*(x-x0)) * np.exp(-(x-x0)**2/(2*sigma**2))
#a =  np.exp(-(x-x0)**2/(2*sigma**2))

# Use dirichlet boundary conditions
a[0] = 0
a[-1] = 0

# could also use np.roll here
ip = np.arange(1,N-1)+1
im = np.arange(1,N-1)-1



#* Initialize plotting variables.
iplot = 1          # Plot counter
aplot = np.copy(a)  # Record the initial state
tplot = np.array([0])       # Record the initial time (t=0)
nplots = 50        # Desired number of plots
plotStep = max(1, np.floor(nStep/nplots)) # Number of steps between plots

a_old = a
a_new = a

#* Loop over desired number of steps.
for iStep in range(1,nStep+1):  ## MAIN LOOP ##
    
    #* Compute new values of wave equation, 
    if iStep == 1:
        a_new[1:N-1] = a[1:N-1] - .5*coeff*(a[ip]-2*a[1:N-1]+a[im])
        a_new[0] = 0
        a_new[N] = 0
        
        a_old = np.copy(a)
        a = np.copy(a_new)
    else:
        a_new[1:N-1] = 2*a[1:N-1] - a_old[1:N-1] + coeff*(a[ip]-2*a[1:N-1]+a[im])
        a_new[0] = 0
        a_new[N] = 0
        
        a_old = np.copy(a)
        a = np.copy(a_new)
  #* Periodically record a(t) for plotting.
    if( (iStep%plotStep) < 1 and iStep >0):  # Every plotStep steps record
        iplot = iplot+1
        aplot = np.vstack((aplot,a))       # Record a(i) for ploting
        tplot = np.append(tplot,tau*iStep)
        #print('%d out of %d steps completed'%(iStep,nStep))

#* Plot the initial and final states.
plt.figure(1)   # Clear figure 1 window and bring forward
plt.plot(x,aplot[0,:],'-',x,a,'--')
plt.legend(['Initial  ','Final'])
plt.xlabel('x')
plt.ylabel('a(x,t)')
plt.grid(True)
# pause(1)    # Pause 1 second between plots

# #* Plot the wave amplitude versus position and time
tt,xx = np.meshgrid(x,tplot)
fig = plt.figure(2)   # Clear figure 2 window and bring forward
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xx, tt, aplot, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
ax.set_xlabel('Time')
ax.set_ylabel('Position')
ax.set_zlabel('Amplitude)')

# mesh(tplot,x,aplot)
# view([-70 50])  # Better view from this angle

plt.show()

# animation
ani = 0
if ani == 1:
    for angle in range(0, 360,4):
        ax.view_init(30, angle)
        plt.draw()
        plt.pause(.001)
