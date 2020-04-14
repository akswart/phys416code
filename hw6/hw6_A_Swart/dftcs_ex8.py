# dftcs - Program to solve the diffusion equation
# using the Forward Time Centered Space (FTCS) scheme.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def runtime(tau):
    plot = False # no plotting
    
    
    #* Initialize parameters (time step, grid spacing, etc.).
    N = 10 # Number of grid points
    L = 1.  # The system extends from x=-L/2 to x=L/2
    h = L/(N-1)  # Grid size
    kappa = 1.   # Diffusion coefficient
    t_sigma = 1/(kappa/h**2*2)
    #print('t_sigma = ',t_sigma)
    #tau = float(input('Enter time step: ')) Get tau from function input
    coeff = kappa*tau/h**2
    """
    if( coeff < 0.5 ):
        print('Solution is expected to be stable')
    else:
        print('WARNING: Solution is expected to be unstable')
    """
    #* Set initial and boundary conditions.
    tt = np.zeros(N)          # Initialize temperature to zero at all points
    tt[round((N-1)/2)] = 1./h     # Initial cond. is delta function in center
    ## The boundary conditions are tt(1) = tt(N) = 0
    
    #* Set up loop and plot variables.
    tt_new = np.zeros(N) # Temp variables for Richardson method
    tt_old = np.zeros(N)
    
    xplot = np.arange(0,N)*h - L/2.   # Record the x scale for plots
    iplot = 1                 # Counter used to count plots
    #nstep = 300               # Maximum number of iterations
    total_time=10
    nstep = int(total_time/tau)
    #print('number of steps =',nstep)
    nplots = 50               # Number of snapshots (plots) to take
    plot_step = nstep/nplots  # Number of time steps between plots
    # record initial time
    ttplot = np.copy(tt)
    tplot = np.array([0.0])
    
    #* Loop over the desired number of time steps.
    
    for istep in range(0,nstep):  ## MAIN LOOP ##
        
        # Start Richardson with FTCS
        if istep == 0:
            #* Compute new temperature using FTCS scheme.
            tt_old = tt.copy()
            tt[1:(N-1)] = tt[1:(N-1)] + coeff*(tt[2:] + tt[0:(N-2)] - 2*tt[1:(N-1)])
            #* Periodically record temperature for plotting.
            if( istep%plot_step < 1 ):   # Every plot_step steps and the first
                ttplot = np.vstack((ttplot,tt))     # record tt(i) for plotting
                tplot = np.append(tplot,istep*tau)      # Record time for plots
                iplot = iplot+1
    
        else: # Continue for n >= 1 with Richardson method
            #* Compute new temperature using Richardson scheme.
            tt_new[1:(N-1)] = tt_old[1:(N-1)] + 2*coeff*(tt[2:] + tt[0:(N-2)] - 2*tt[1:(N-1)])
            tt_old = tt.copy()
            tt = tt_new.copy()
            
            #* Periodically record temperature for plotting.
            if( istep%plot_step < 1 ):   # Every plot_step steps and the first
                ttplot = np.vstack((ttplot,tt))     # record tt(i) for plotting
                tplot = np.append(tplot,istep*tau)      # Record time for plots
                iplot = iplot+1
        if max(tt) >= 10./h:
             # If we blowup past 10 times initial delta func value then 
             # its clear we won't converge so return last value array and time to blowup
            return tt,istep*tau,istep # If we blowup past initial delta func value
        
    if plot:
        #* Plot temperature versus x and t as wire-mesh and contour plots.
        tt,xx = np.meshgrid(xplot,tplot)
        fig = plt.figure(1)
        plt.clf()
        ax = fig.gca(projection='3d')
        # mesh(tplot,xplot,ttplot)  # Wire-mesh surface plot
        surf = ax.plot_surface(xx, tt, ttplot, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
        ax.set_xlabel('Time')
        ax.set_ylabel('x')
        ax.set_zlabel('T(x,t)')
        # ax.set_title('Diffusion of a delta spike')
        # pause(1)
        
        plt.figure(2)
        plt.clf()
        contourLevels = np.arange(0,10,0.5) #  contourLabels = 0:5
        plt.contour(xx,tt,ttplot,contourLevels)  # Contour plot
        # clabel(cs,contourLabels)  # Add labels to selected contour levels
        plt.xlabel('Time')
        plt.ylabel('x')
        plt.title('Temperature contour plot')
        plt.show()
        
    return tt_old,-1,-1

if __name__ == '__main__':
    N = 10 # Number of grid points
    L = 1.  # The system extends from x=-L/2 to x=L/2
    h = L/(N-1)  # Grid size
    kappa = 1.   # Diffusion coefficient
    t_max = 1/(kappa/h**2*2)
    t_scale = 17
    
    step_list = []
    tt_val = []
    t_blowup = []
    tau_list = [t_max/(2**i) for i in range(1,t_scale)]
    for tau in tau_list:
        tt,time,istep = runtime(tau)
        
        t_blowup.append(time)
        tt_val.append(tt)
        step_list.append(istep)
    
    plt.figure(1)
    plt.semilogx(tau_list,t_blowup)
    plt.xlabel('Tau')
    plt.ylabel('Time to blow up')
    #print(t_blowup)
    #print(tau_list)
    #print(step_list)
    plt.figure(2)
    plt.loglog(tau_list,step_list)
    plt.xlabel('Tau')
    plt.ylabel('Steps to blow up')

