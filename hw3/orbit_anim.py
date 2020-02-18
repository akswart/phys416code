import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# python 3 version 2/15
def rk4(x,t,tau,derivsRK):
#%  Runge-Kutta integrator (4th order)
#% Input arguments -
#%   x = current value of dependent variable
#%   t = independent variable (usually time)
#%   tau = step size (usually timestep)
#%   derivsRK = right hand side of the ODE; derivsRK is the
#%             name of the function which returns dx/dt
#%             Calling format derivsRK(x,t).
#% Output arguments -
#%   xout = new value of x after a step of size tau
    half_tau = 0.5*tau
    F1 = derivsRK(x,t)
    t_half = t + half_tau
    xtemp = x + half_tau*F1
    F2 = derivsRK(xtemp,t_half)
    xtemp = x + half_tau*F2
    F3 = derivsRK(xtemp,t_half)
    t_full = t + tau
    xtemp = x + tau*F3
    F4 = derivsRK(xtemp,t_full)
    xout = x + tau/6.*(F1 + F4 + 2.*(F2+F3))
    return xout
def  rka(x,t,tau,err,derivsRK):

#% Adaptive Runge-Kutta routine
#% Inputs
#%   x          Current value of the dependent variable
#%   t          Independent variable (usually time)
#%   tau        Step size (usually time step)
#%   err        Desired fractional local truncation error
#%   derivsRK   Right hand side of the ODE; derivsRK is the
#%              name of the function which returns dx/dt
#%              Calling format derivsRK(x,t).
#% Outputs
#%   xSmall     New value of the dependent variable
#%   t          New value of the independent variable
#%   tau        Suggested step size for next call to rka

#%* Set initial variables
    tSave = t;  xSave = x    # Save initial values
    safe1 = .9;  safe2 = 4.  # Safety factors
    eps = np.spacing(1) # smallest value

#%* Loop over maximum number of attempts to satisfy error bound
    maxTry = 100

    for iTry in range(1,maxTry):

#%* Take the two small time steps
        half_tau = 0.5 * tau
        xTemp = rk4(xSave,tSave,half_tau,derivsRK)
        t = tSave + half_tau
        xSmall = rk4(xTemp,t,half_tau,derivsRK)

  #%* Take the single big time step
        t = tSave + tau
        xBig = rk4(xSave,tSave,tau,derivsRK)

  #%* Compute the estimated truncation error
        scale = err * (np.abs(xSmall) + np.abs(xBig))/2.
        xDiff = xSmall - xBig
        errorRatio = np.max( [np.abs(xDiff)/(scale + eps)] )

        #print safe1,tau,errorRatio

  #%* Estimate news tau value (including safety factors)
        tau_old = tau

        tau = safe1*tau_old*errorRatio**(-0.20)
        tau = np.max([tau,tau_old/safe2])
        tau = np.min([tau,safe2*tau_old])

  #%* If error is acceptable, return computed values
        if errorRatio < 1 :
            xSmall = xSmall #% +  (xDiff)/15
            return xSmall, t, tau

#%* Issue error message if error bound never satisfied
    print ('ERROR: Adaptive Runge-Kutta routine failed')
    return
def gravrk(s,t):
#%  Returns right-hand side of Kepler ODE; used by Runge-Kutta routines
#%  Inputs
#%    s      State vector [r(1) r(2) v(1) v(2)]
#%    t      Time (not used)
#%  Output
#%    deriv  Derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]

#%* Compute acceleration
    r = np.array([s[0], s[1]])  # Unravel the vector s into position and velocity
    v = np.array([s[2] ,s[3]])
    accel = -GM*r/np.linalg.norm(r)**3  # Gravitational acceleration

#%* Return derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]
    derivs = np.array([v[0], v[1], accel[0], accel[1]])
    return derivs

# orbit - Program to compute the orbit of a comet.
#clear all;  help orbit;  % Clear memory and print header

# Set initial position and velocity of the comet.
r0 = float(input("Enter initial radial distance (AU): "))
# v0 = float(input("Enter initial tangential velocity (AU/yr): "))
# r0 = 1.0; v0 = np.pi;
#jag20 2.5.15 modify input to allow 'pi'
vstr = str(input("Enter initial tangential velocity (AU/yr) as decimal or decimal multiple of Pi (e.g.,2*pi): "))
vinp = vstr.split('*')
if (vinp[-1].lower() == 'pi'):
  v0 = float(vinp[0])*np.pi
else:
  v0 = float(vinp[0])

r = np.array([r0, 0.]);  v = np.array([0., v0])
state = np.array([ r[0], r[1], v[0], v[1] ])   # Used by R-K routines

#Set physical parameters (mass, G*M)
GM = 4*np.pi**2      # Grav. const. * Mass of Sun (au^3/yr^2)
mass = 1.        # Mass of comet
adaptErr = 1.e-4 # Error parameter used by adaptive Runge-Kutta
time = 0.0

energy_factor = 2/r0 - v0**2/GM;
if(energy_factor <0):
    print('warning, orbit not bound')
else:
    period = 1./(2./r0 - v0**2/GM)**1.5;
    print('Period or the orbit = %f years' %period)

#%* Loop over desired number of steps using specified
#%  numerical method.
nStep = int(input("Enter number of steps: "))
tau =  float(input("Enter time step (yr): "))
# fixme: should put in a statement here to check that an acceptable number was chosen
NumericalMethod = int(input("Choose a number for a numerical method:  \
1-Euler, 2-Euler-Cromer, 3-Runge-Kutta 4-Adaptive R-K: "))

for istep in range(0,nStep):

  #%* Record position and energy for plotting.
  # Initially set the arrays for the first step
  if istep == 0:
      rplot = np.linalg.norm(r)
      thplot = np.arctan2(r[1],r[0])
      tplot = time
      kinetic = .5*mass*np.linalg.norm(v)**2
      potential= - GM*mass/np.linalg.norm(r)

  rplot = np.append(rplot,np.linalg.norm(r))           #Record position for polar plot
  thplot = np.append(thplot,np.arctan2(r[1],r[0]))
  tplot = np.append(tplot,time)
  kinetic = np.append(kinetic,0.5*mass*np.linalg.norm(v)**2)   # Record energies
  potential= np.append(potential,- GM*mass/np.linalg.norm(r))

  #%* Calculate new position and velocity using desired method.
  if NumericalMethod == 1 :
    accel = -GM*r/np.linalg.norm(r)**3
    r = r + tau*v             # Euler step
    v = v + tau*accel
    time = time + tau
  elif NumericalMethod == 2:
    accel = -GM*r/np.linalg.norm(r)**3
    v = v + tau*accel
    r = r + tau*v             # Euler-Cromer step
    time = time + tau
  elif NumericalMethod == 3:
    state = rk4(state,time,tau,gravrk)
    r = [state[0], state[1]]   # 4th order Runge-Kutta
    v = [state[2], state[3]]
    time = time + tau
  else:
    [state, time, tau] = rka(state,time,tau,adaptErr,gravrk)
    r = [state[0], state[1]]   # Adaptive Runge-Kutta
    v = [state[2], state[3]]

#%* Graph the trajectory of the comet.
plt.figure(1); plt.clf()  #Clear figure 1 window and bring forward
plt.polar(thplot,rplot,'-')  # Use polar plot for graphing orbit
plt.xlabel('Distance (AU)')
plt.grid(True)

#%* Graph the energy of the comet versus time.
plt.figure(2); plt.clf()   # Clear figure 2 window and bring forward
totalE = kinetic + potential   # Total energy
plt.plot(tplot,kinetic,'-.',tplot,potential,'--',tplot,totalE,'-')
#plt.legend('Kinetic','Potential','Total')
plt.xlabel('Time (yr)'); plt.ylabel('Energy (M AU^2/yr^2)')
plt.grid(True)
plt.show()

#plt.ion() # turn on interactive

"""
# animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)
"""

xearth = rplot*np.cos(thplot)
yearth = rplot*np.sin(thplot)
xsun = 0.0; ysun = 0.0;

fig = plt.figure(3);
plt.axis('equal')
plt.axis([-1.5,1.5,-1.5,1.5])

"""
def animate(i):
    plt.title('time =%0.2f years'%tplot[i])
    plt.plot(xearth[i],yearth[i],'bo')
    plt.plot(xsun,ysun,'y*')

#plt.title('time =%0.2f years'%tplot[i])
ani = animation.FuncAnimation(fig, animate, frames=len(xearth), repeat=True)
#ani.save('HeroinOverdosesJumpy.mp4', writer=writer)
plt.show()

"""
for i in range(0,len(xearth)):
    plt.clf() # clear
    plt.plot(xearth[0:i],yearth[0:i],'b-')
    plt.plot(xearth[i],yearth[i],'bo')
    plt.plot(xsun,ysun,'y*')
    plt.axis('equal')
    plt.axis([-1.5,1.5,-1.5,1.5])
    plt.title('time =%0.2f years'%tplot[i])
    plt.draw()
    plt.pause(tau)
    plt.show()
#"""