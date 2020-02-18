#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 17:01:26 2020

@author: aswart
"""

# python 3 version 2/15
import numpy as np
import matplotlib.pyplot as plt

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
def rka(x,t,tau,err,derivsRK):

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
          # xSmall = xSmall #% +  (xDiff)/15
          #   xSmall = (16.*xSmall - xBig)/15. # correction
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

    GM = 4*np.pi**2
    
#%* Compute acceleration
    r = np.array([s[0], s[1]])  # Unravel the vector s into position and velocity
    v = np.array([s[2] ,s[3]])
    accel = -GM*r/np.linalg.norm(r)**3  # Gravitational acceleration

#%* Return derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]
    derivs = np.array([v[0], v[1], accel[0], accel[1]])
    return derivs

# orbit - Program to compute the orbit of a comet.
"""
Planning:
    Input:
    Takes list of objects
    Each object has a dict that states:
        radius (r0)
        velocity (v0)
        mass (mass)
    tau and nSteps is directly passed in
    
    Output:
    Iterate over array and build inital value arrays for each planet
    each planet is a dict in a list
    dict contains variables with standard names as from orbit.py as keys
    values for each key are lists
    
"""
def orbit(tau, nStep, input_list = [],calc_info = False, plot_momentum = False,
           plot_traj = True, plot_energy = False):        

    if input_list:
        for planet in input_list:
            print(planet)
    else:
        r0_e = 1
        v0_e = 2*np.pi
        
        r0_j = 4
        v0_j = 2*np.pi

    
    r_e = np.array([r0_e, 0.])
    v_e = np.array([0., v0_e])
    state_e = np.array([ r_e[0], r_e[1], v_e[0], v_e[1] ])   # Used by R-K routines
    
    r_j = np.array([r0_j, 0.])
    v_j = np.array([0., v0_j])
    state_j = np.array([ r_j[0], r_j[1], v_j[0], v_j[1] ]) 
    
    #Set physical parameters (mass, G*M)
    GM = 4*np.pi**2      # Grav. const. * Mass of Sun (au^3/yr^2)
    mass_e = 1.        # Mass of comet
    mass_j = 1.
    adaptErr = 1.e-4 # Error parameter used by adaptive Runge-Kutta
    time = 0.0
    
    #%* Loop over desired number of steps using specified
    #%  numerical method.
    for istep in range(0,nStep):
    
      #%* Record position and energy for plotting.
      # Initially set the arrays for the first step
      if istep == 0:
          tplot = time

          rplot_e = np.linalg.norm(r_e)
          thplot_e = np.arctan2(r_e[1],r_e[0])
          kinetic_e = .5*mass_e*np.linalg.norm(v_e)**2
          potential_e = - GM*mass_e/np.linalg.norm(r_e)
          momentum_e = [np.linalg.norm(np.cross(r_e,mass_e*v_e))]
          
          rplot_j = np.linalg.norm(r_j)
          thplot_j = np.arctan2(r_j[1],r_j[0])
          kinetic_j = .5*mass_j*np.linalg.norm(v_j)**2
          potential_j = - GM*mass_j/np.linalg.norm(r_j)
          momentum_j = [np.linalg.norm(np.cross(r_j,mass_j*v_j))]
      else:    
          tplot = np.append(tplot,time)
          
          rplot_e = np.append(rplot_e,np.linalg.norm(r_e))           #Record position for polar plot
          thplot_e = np.append(thplot_e,np.arctan2(r_e[1],r_e[0]))
          kinetic_e = np.append(kinetic_e,0.5*mass_e*np.linalg.norm(v_e)**2)   # Record energies
          potential_e = np.append(potential_e,- GM*mass_e/np.linalg.norm(r_e))
          momentum_e.append(np.linalg.norm(np.cross(r_e, mass_e*v_e)))
          
          rplot_j = np.append(rplot_j,np.linalg.norm(r_j))           #Record position for polar plot
          thplot_j = np.append(thplot_j,np.arctan2(r_j[1],r_j[0]))
          kinetic = np.append(kinetic_j,0.5*mass_j*np.linalg.norm(v_j)**2)   # Record energies
          potential= np.append(potential_j,- GM*mass_j/np.linalg.norm(r_j))
          momentum_j.append(np.linalg.norm(np.cross(r_j, mass_j*v_j)))

      #%* Calculate new position and velocity using Adaptive RK4
        [state, time, tau] = rka(state,time,tau,adaptErr,gravrk)
        r = np.array([state[0], state[1]])   # Adaptive Runge-Kutta
        v = np.array([state[2], state[3]])
      # If sometime after first step and theta goes from neg to pos, then we know we've completed an orbit
 #     if istep != 0 and (thplot[-2:]*[-1,1] > 0).all():
 #         break
          
    #print(thplot[-5:])
    
    if plot_traj:
        #%* Graph the trajectory of the comet.
        plt.figure(1); plt.clf()  #Clear figure 1 window and bring forward
        plt.polar(thplot,rplot,'-')  # Use polar plot for graphing orbit
        plt.xlabel('Distance (AU)')
        plt.grid(True)
    
    if plot_energy:
        #%* Graph the energy of the comet versus time.
        plt.figure(2); plt.clf()   # Clear figure 2 window and bring forward
        totalE = kinetic + potential   # Total energy
        plt.plot(tplot,kinetic,'-.',tplot,potential,'--',tplot,totalE,'-')
        #plt.legend('Kinetic','Potential','Total')
        plt.xlabel('Time (yr)'); plt.ylabel('Energy (M AU^2/yr^2)')
        plt.grid(True)
        plt.show()
        
    if plot_momentum:
        # Plots angular momentum as a function of time
        plt.figure(3)
        plt.plot(tplot,momentum)
   
    

    return rplot, thplot

if __name__ == "__main__":
    # non-elliptical
    input_dict = {
        'r0': 1,
        'v0': 1*np.pi,
        'nStep': 1000,
        'tau': .01,
        }
    
    
    
    rplot, thplot  = orbit(.01,1000, input_dict)
