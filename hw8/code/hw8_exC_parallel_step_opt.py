import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from timeit import default_timer as timer

"""
RUN AT YOUR OWN RISK, THIS IS MULTITHREADED BENCHMARKING DONE WITH 
THE CODE FROM THE LAST PROBLEM. IT WILL USE ALL CPU RESOURCES AND CAN
UPSET ANTIVIRUS PROGRAMS, (THOUGH I MANAGED TO MAKE THAT ISSUE GO AWAY FOR ME)

IF YOU WANT TO RUN IT, I SUGGEST THE FOLLOWING SETTINGS AT THE BOTTOM:
    beg_scale = .5
    end_scale = 1
    num_point = 5
IN THE BOTTOM IF-STATEMENT. THESE ARE NOT THE DEFAULT VALUES.
"""


def corona(N_scale):
    start = timer()
    #* Initialize parameters (system size, grid spacing, etc.)
    method = 3
    rN = int(61*N_scale)
    pN = int(91*N_scale)
    animate=0
    graph =0
    n = 3
    rL = 2          # System size (length)
    pL = 2*np.pi
    rh = rL/(rN-1)    # Grid spacing
    ph = pL/(pN-1)
    r = np.arange(0,rN)*rh+1  # r coordinate
    phi = np.arange(0,pN)*ph # phi coordinate
    x = np.array([i*np.cos(j) for i in r for j in phi])
    y = np.array([i*np.sin(j) for i in r for j in phi])
    R,P = np.meshgrid(r,phi) # for plotting, note the reversal in x and y
    X, Y = R*np.cos(P), R*np.sin(P)
    
    plot_interval = 50 # interval to plot animation, setting it smaller slows the program down alot
    
    #* Select over-relaxation factor (SOR only)
    if( method == 3 ):
        rad = .5*(np.cos(np.pi/rN) + np.cos(np.pi/pN))
        omegaOpt = 2/(1+np.sqrt(1-rad**2))  # Theoretical optimum
        #print('Theoretical optimum omega = ',omegaOpt)
        omega = omegaOpt#float(input('Enter desired omega: '))
    
    #* Set initial guess as first term in separation of variables soln.
    A0 = 1     # Potential at r=1
    # phi = phi0 * 4/(np.pi*np.sinh(np.pi)) * np.outer(np.sin(np.pi*x/L),np.sinh(np.pi*y/L))
    A=np.zeros((rN,pN)) # try this to see it evolve better
    
    #* Set boundary conditions
    # first index is the radius and second index is phi (rows,cols)
    for i in range(pN): # Apply inner boundary conditions 
        A[0,i] = np.cos(n*phi[i])
    A[-1,:] = A[-2,:] # Apply outer boundary condition
    A[:,0] = A[:,-1] # Periodic Boundary conditions about phi = 0
    
    
    #print('Potential is zero on all other boundaries')
    
    #plt.ion()
    
    #* Loop until desired fractional change per iteration is obtained
    # start_time=cputime     # Reset the cputime counter
    newphi = np.copy(phi)           # Copy of the solution (used only by Jacobi)
    iterMax = pN**2          # Set max to avoid excessively long runs
    changeDesired = 1.0e-4   # Stop when the change is given fraction
    #print('Desired fractional change = ',changeDesired)
    change = np.array([])
    
    for iterk in range(0,iterMax):
        changeSum = 0.0
      
        ## SOR method ##
        for i in range(1,rN-1):        # Loop over interior points only
            for j in range(1,pN-1):
                newA = 0.25*omega*(A[i+1,j]+A[i-1,j]+ A[i,j-1]+A[i,j+1])  +  (1-omega)*A[i,j]
                changeSum = changeSum + abs(1-A[i,j]/newA)
                A[i,j] = newA
        
        
        # Update boundary conditions
        for i in range(pN): # Apply inner boundary conditions 
            A[0,i] = np.cos(n*phi[i])
        A[-1,:] = A[-2,:] # Apply outer boundary condition
        A[:,0] = A[:,-1] # Periodic Boundary conditions about phi = 0
        
        #* Check if fractional change is small enough to halt the iteration
        change = np.append(change,changeSum/(pN-2)**2)
        #if( iterk%10 < 1 ):
            #print('After %d iterations, fractional change = %f'%( iter,change[-1]))
            
    
        if( change[-1] < changeDesired ):
          #print('Desired accuracy achieved after %d iterations'%iter)
          #print('Breaking out of main loop')
          break
    # animate
        if(animate ==1 and iterk%plot_interval<1):
            fig = plt.figure(2)   # Clear figure 2 window and bring forward
            plt.clf()
            ax = fig.gca(projection='3d')
            surf = ax.plot_surface(X, Y, A.T, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('potential after '+str(iterk)+' iterations')
            plt.draw()
            plt.show()
            plt.pause(0.1)
    
    # total_time = cputime - start_time # get the total cpu time
    
    #* Plot final estimate of potential as contour and surface plots
    
    #plt.ioff()
    if graph ==1:
        plt.figure(1);plt.clf()
        contourLevels = np.arange(0,1,0.1) #
        plt.contour(X,Y,A.T,contourLevels)  # Contour plot
        # clabel(cs,contourLabels)  # Add labels to selected contour levels
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(r'$\Phi(x,y)$')
        
        fig = plt.figure(2)   # Clear figure 2 window and bring forward
        plt.clf()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X, Y, A.T, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('potential after '+str(iterk)+' iterations')
        
        #* Plot the fractional change versus iteration
        plt.figure(3); 
        plt.clf()
        plt.semilogy(change)
        plt.xlabel('Iteration')
        plt.ylabel('Fractional change')
        if method==1:
            title=' Jacobi method'
        elif(method==2):
            title =' Gauss-Seidel method'
        elif(method==3):
            title=' SOR method, $\Omega$ ='+str(omega)
        
        plt.title(r'Iterations ='+str(iterk)+title+' pN='+str(pN))
        plt.grid(True)
        
        fig = plt.figure(4)   # Clear figure 2 window and bring forward
        plt.clf()
        ax = fig.gca(projection='3d')
        A_analytic = np.zeros((rN,pN))
        for i in range(rN):
            for j in range(pN):
                A_analytic[i,j] = ( ( (3**(2*n)/r[i]**2) + r[i]**n)/(3**(2*n)   ) )*np.cos(n*phi[j])
        surf = ax.plot_surface(X, Y, A_analytic.T, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('potential from analytic soln')
        plt.show()
    stop = timer()
    time = stop-start
    return (iterk,time)
    
"""
Lets benchmark convergence speed for a bunch of different grid sizes

"""

if __name__ == "__main__":
    import multiprocessing
    
    beg_scale = .5
    end_scale = 1
    num_point = 5
    
    times = [] # Append tuple of (scale, time, interations) 
    
    scales = np.linspace(beg_scale,end_scale,num_point)
    
    # Start parallel pool
    p = multiprocessing.Pool()
    result = p.map(corona,scales)
    p.close()
    p.join()
    for i in range(len(result)):
        times.append( (scales[i], result[i][1], result[i][0]) )
    
    import pandas as pd
    df = pd.DataFrame(times)
    df.columns = ["Scale","Time","Iterations"]
    df.to_csv("outputfile.csv",index=False)    








