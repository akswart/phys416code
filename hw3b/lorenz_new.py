# python 3 version
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

def rk4(x,t,tau,derivsRK):
    """
    ##  Runge-Kutta integrator (4th order)
    ## Input arguments -
    ##   x = current value of dependent variable
    ##   t = independent variable (usually time)
    ##   tau = step size (usually timestep)
    ##   derivsRK = right hand side of the ODE; derivsRK is the
    ##             name of the function which returns dx/dt
    ##             Calling format derivsRK(x,t).
    ## Output arguments -
    ##   xout = new value of x after a step of size tau
    """
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
    """
    ## Adaptive Runge-Kutta routine
    ## Inputs
    ##   x          Current value of the dependent variable
    ##   t          Independent variable (usually time)
    ##   tau        Step size (usually time step)
    ##   err        Desired fractional local truncation error
    ##   derivsRK   Right hand side of the ODE; derivsRK is the
    ##              name of the function which returns dx/dt
    ##              Calling format derivsRK(x,t).
    ## Outputs
    ##   xSmall     New value of the dependent variable
    ##   t          New value of the independent variable
    ##   tau        Suggested step size for next call to rka
    """
    # Set initial variables
    tSave = t;  xSave = x    # Save initial values
    safe1 = .9;  safe2 = 4.  # Safety factors
    eps = np.spacing(1) # smallest value

    # Loop over maximum number of attempts to satisfy error bound
    maxTry = 100

    for iTry in range(1,maxTry):
	
        # Take the two small time steps
        half_tau = 0.5 * tau
        xTemp = rk4(xSave,tSave,half_tau,derivsRK)
        t = tSave + half_tau
        xSmall = rk4(xTemp,t,half_tau,derivsRK)
          
        # Take the single big time step
        t = tSave + tau
        xBig = rk4(xSave,tSave,tau,derivsRK)
          
        # Compute the estimated truncation error
        scale = err * (np.abs(xSmall) + np.abs(xBig))/2.
        xDiff = xSmall - xBig
        errorRatio = np.max( [np.abs(xDiff)/(scale + eps)] )
        
        #print safe1,tau,errorRatio
          
        # Estimate news tau value (including safety factors)
        tau_old = tau
        
        tau = safe1*tau_old*errorRatio**(-0.20)
        tau = np.max([tau,tau_old/safe2])
        tau = np.min([tau,safe2*tau_old])
  
    # If error is acceptable, return computed values
    if errorRatio < 1 : 
        xSmall = xSmall 
        return xSmall, t, tau  
    # Issue error message if error bound never satisfied
    print ('ERROR: Adaptive Runge-Kutta routine failed')
    return



def lorzrk(s,t):
    """
    #  Returns right-hand side of Lorenz model ODEs
    #  Inputs
    #    s      State vector [x y z]
    #    t      Time (not used)
    #    param  Parameters [r sigma b]
    #  Output
    #    deriv  Derivatives [dx/dt dy/dt dz/dt]
    """
    #* For clarity, unravel input vectors
    x = s[0]; y = s[1]; z = s[2]
    #* Return the derivatives [dx/dt dy/dt dz/dt]
    deriv = np.zeros(3)
    deriv[0] = sigma*(y-x)
    deriv[1] = r*x - y - x*z
    deriv[2] = x*y - b*z
    return deriv

# Set initial state x,y,z and parameters r,sigma,b

sxin, syin, szin = input("Enter the initial position x, y, z: ").split(',')
state = np.zeros(3)
state[0] = float(sxin); state[1]  = float(syin); state[2]  = float(szin)

r = float(input('Enter the parameter r: '))

sigma = 10   # Parameter sigma
b = 8./3.     # Parameter b
param = np.array([r, sigma, b])  # Vector of parameters passed to rka
tau = 1       # Initial guess for the timestep
err = 1.e-3   # Error tolerance

# Loop over the desired number of steps
time = 0;
nstep = int(input('Enter number of steps: '))
# initialize arrays
tplot=np.array([]); tauplot=np.array([])
xplot=np.array([]); yplot=np.array([]); zplot=np.array([])

for istep in range(0,nstep):
    # Record values for plotting
    x = state[0]; y = state[1]; z = state[2];
    tplot = np.append(tplot,time)
    tauplot = np.append(tauplot,tau)
    xplot = np.append(xplot,x)
    yplot = np.append(yplot,y)
    zplot = np.append(zplot,z)
    if( istep%50 ==0 ):
        print('Finished %d steps out of %d '%(istep,nstep))

    #* Find new state using adaptive Runge-Kutta
    [state, time, tau] = rka(state,time,tau,err,lorzrk)

# Print max and min time step returned by rka
print('Adaptive time step: Max = %f,  Min = %f '%(max(tauplot[1:]), min(tauplot[1:])));
    
#* Graph the time series x(t)
plt.figure(1); plt.clf();  # Clear figure 1 window and bring forward
plt.plot(tplot,xplot,'-')
plt.xlabel('Time');  plt.ylabel('x(t)')
plt.title('Lorenz model time series')
# plt.show()

# #* Graph the x,y,z phase space trajectory

fig=plt.figure(2)
ax=p3.Axes3D(fig)
ax.plot3D(xplot,yplot,zplot)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.grid(True)
# title('Lorenz model phase space')
ax.set_aspect('equal')
plt.show()




