# traffic - Program to solve the generalized Burger  
# equation for the traffic at a stop light problem
# clear all  help traffic # Clear memory and print header
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm # colormap
from mpl_toolkits.mplot3d import Axes3D

#* Select numerical parameters (time step, grid spacing, etc.).
#method = int(input('Choose a numerical method: (1)FTCS, (2)Lax, (3)Lax-Wendroff: '))
method = 3
#N = int(input('Enter number of grid points: '))
N = 80
L = 400.     # System size (meters)
h = L/N      # Grid spacing for periodic boundary conditions
v_max = 25.   # Maximum car speed (m/s)
print('Suggested timestep is ',h/v_max)
#tau = float(input('Enter time step (tau): '))
tau = h/v_max
print('Last car starts moving after %f steps'%((L/4)/(v_max*tau)))
#nstep = int(input('Enter number of steps: '))
nstep = 30
coeff = tau/(2*h)        # Coefficient used by all schemes
coefflw = tau**2/(2*h**2)  # Coefficient used by Lax-Wendroff

#* Set initial and boundary conditions
rho_max = 1.0                  # Maximum density
Flow_max = 0.25*rho_max*v_max  # Maximum Flow

#x = np.arange()


""" # Old init conds
# Initial condition is a square pulse from x = -L/4 to x = 0
rho = np.zeros(N)

for i in range(int(N/4),int(N/2-1)):
  rho[i] = rho_max    # Max density in the square pulse

rho[int(N/2-1)] = rho_max/2  # Try running without this line
"""
# Use periodic boundary conditions
ip = np.arange(0,N)+1
ip[N-1] = 0 # ip = i+1 with periodic b.c.
im = np.arange(0,N)-1
im[0] = N-1 # im = i-1 with periodic b.c.

#* Initialize plotting variables.
iplot = 1
xplot = (np.arange(0,N)+1./2.)*h - L/2  # Coordinates of grid points

# New initial conditions
rho = (rho_max/2)*(1+np.cos(4*np.pi*xplot/L))

rplot = np.copy(rho)       # Record the initial state
tplot = np.array([0.0])
plt.ion()
plt.figure(1)
plt.clf()  # Clear figure 1 window and bring forward

#* Loop over desired number of steps.
for istep in range(0,nstep):
	
  #* Compute the flow = (Density)*(Velocity)
  Flow = rho * (v_max*(1 - rho/rho_max))

  #* Compute new values of density using FTCS, 
  #  Lax or Lax-Wendroff method.
  if( method == 1 ) :     ### FTCS method ###
    rho[0:N] = rho[0:N] - coeff*(Flow[ip]-Flow[im])
  elif( method == 2 ):  ### Lax method ###
    rho[0:N] = 0.5*(rho[ip]+rho[im]) - coeff*(Flow[ip]-Flow[im])
  else:                   ### Lax-Wendroff method ###
    cp = v_max*(1. - (rho[ip]+rho[0:N])/rho_max)
    cm = v_max*(1. - (rho[0:N]+rho[im])/rho_max)
    rho[0:N] = rho[0:N] - coeff*(Flow[ip]-Flow[im]) + \
               coefflw*(cp*(Flow[ip]-Flow[0:N]) - cm*(Flow[0:N]-Flow[im]))

  #* Record density for plotting.
  iplot = iplot+1
  rplot = np.vstack((rplot,rho))
  tplot = np.append(tplot,tau*(istep+1))
  
  #* Display snap-shot of density versus position
  plt.clf()
  plt.plot(xplot,rho,'-',xplot,Flow/Flow_max,'--')
  plt.xlabel('x')
  plt.ylabel('Density and Flow')
  plt.title(' time ='+"{:12.2f}".format(tplot[-1])+'sec')
  plt.legend([r'$\rho(x,t)$','F(x,t)'])
  plt.axis([-L/2, L/2, -0.1, 1.1])
  plt.draw()
  plt.pause(0.1)


plt.ioff()
#* Graph density versus position and time as wire-mesh plot
tt,xx = np.meshgrid(xplot,tplot)
fig = plt.figure(2)   # Clear figure 2 window and bring forward
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xx, tt, rplot, rstride=1, cstride=1) #, cmap=cm.jet,linewidth=0, antialiased=False)
ax.set_xlabel('Time')
ax.set_ylabel('Position')
ax.set_zlabel('Density versus position and time)')

plt.figure(3)
contourLevels = np.arange(0,1,0.1) #
plt.contour(tt,xx,rplot,contourLevels)  # Contour plot
# clabel(cs,contourLabels)  # Add labels to selected contour levels
plt.xlabel('Time')
plt.ylabel('x')
plt.title('Density contours')
plt.show()


