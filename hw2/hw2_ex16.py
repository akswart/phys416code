# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 10:04:27 2020

@author: akswa
"""

# pendul - Program to compute the motion of a simple pendulum
# using the Euler or Verlet method

# Select the numerical method to use: Euler or Verlet


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk








def period_pend(theta0,g_over_L):
# function to return the exact period for a pendulum of length L
# usage: period = exact_period(theta0,g_over_L)
# where: theta0 = inital angle in degrees
#        g_over_L = ratio g to the length of the pendulum
#        note  -earlier version has a bug as it x sqrt(g/l) not divided 9/11

# note the squaring of the argument in the elliptic function
# matlab uses a different normalization than the book

    period = 4/np.sqrt(g_over_L)*ellipk((np.sin(theta0*np.pi/180./2.))**2)
    return period

NumericalMethod = int(input("Choose a numerical method: Euler(1) Verlet(2): "))


def pend(theta0,tau,nstep,NumericalMethod,plotting = False):
    """
    Parameters
    ----------
    theta0 : Float
        Intitial Angle of pendulum.
    tau : Float
        Timestep size (seconds).
    nstep : Int
        Number of timesteps.
    nstep : string
        Type of numerical method.
        Options: "euler", "verlet", "eulercromer", "leapfrog", "midpoint" 
    Returns
    -------
    None.
    
    """	   
    
    # Set initial position and velocity of pendulum
    theta = theta0*np.pi/180  # Convert angle to radians
    omega = 0               # Set the initial velocity
    
    # Set the physical constants and other variables
    g_over_L = 1            # The constant g/L
    time = 0                # Initial time
    irev = 0                # Used to count number of reversals
    
    # Take one backward step to start Verlet
    accel = -g_over_L*np.sin(theta)    # Gravitational acceleration
    theta_old = theta - omega*tau + 0.5*tau**2*accel    
    
    # Loop over desired number of steps with given time step
    #    and numerical method
    #nstep = int(input("Enter number of time steps: "))
    # initialize arrays
    t_plot=np.array([])
    th_plot=np.array([])
    period=np.array([])
    
    for istep in range(0,nstep):  
    
        # Record angle and time for plotting
        t_plot = np.append(t_plot,time)            
        th_plot = np.append(th_plot,theta*180/np.pi)   # Convert angle to degrees
        time = time + tau
      
        # Compute new position and velocity using Euler or Verlet method
        accel = -g_over_L*np.sin(theta);    # Gravitational acceleration
        if NumericalMethod == "euler":
            if istep == 0: print("using euler method")
            theta_old = theta               # Save previous angle
            theta = theta + tau*omega       # Euler method
            omega = omega + tau*accel 
        else: 
            if istep == 0: print("using verlet method")
            theta_new = 2*theta - theta_old + tau**2*accel
            theta_old = theta			   # Verlet method
            theta = theta_new  
      
        # Test if the pendulum has passed through theta = 0;
        #    if yes, use time to estimate period
        if theta*theta_old < 0: # Test position for sign change
            print("Turning point at time t= %f" %time) ;
            if irev == 0:           # If this is the first change,
                time_old = time    # just record the time
            else:
                period = np.append(period,2*(time - time_old))
                time_old = time
            irev = irev + 1       # Increment the number of reversals
    
    if irev > 1:
    # Estimate period of oscillation, including error bar
        AvePeriod = np.mean(period)
        ErrorBar = np.std(period)/np.sqrt(irev)
        print("Average period = %g +/- %g" %(AvePeriod,ErrorBar))
    else:
        print('Pendulum program could not complete a period, time =%g'%time)
    
    print("Exact period = %g" %period_pend(theta0,g_over_L))
    
    # Graph the oscillations as theta versus time
    if plotting:
        plt.figure(1)         # Clear and forward figure window
        plt.plot(t_plot,th_plot,'+')
        plt.xlabel('Time')
        plt.ylabel(r'$\theta$ (degrees)') # the 'r' means raw strings for latex
        plt.grid()
        plt.show()
    
    return t_plot,th_plot,period


if __name__ == "__main__":
    a,b,c = pend(10,.01,1000,"euler",True)
    
