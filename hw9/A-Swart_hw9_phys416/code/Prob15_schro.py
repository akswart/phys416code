import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# from mpl_toolkits.mplot3d import Axes3D
#  schro - Program to solve the Schrodinger equation
#  for a free particle using the Crank-Nicolson scheme
# clear all  help schro   # Clear memory and print header

def schro(U):

    #* Initialize parameters (grid spacing, time step, etc.)
    i_imag = 1j #np.sqrt(-1)    # Imaginary i
    #N = int(input('Enter number of grid points: '))
    N = 2001
    L = 600              # System extends from -L/2 to L/2
    h = L/(N-1)          # Grid size
    x = (np.arange(0,N))*h - L/2  # Coordinates of grid points
    h_bar = 1
    mass = 1 # Natural units
    
    #tau = float(input('Enter time step: '))
    tau = .05
    
    #* Initialize the wavefunction 
    x0 = 0          # Location of the center of the wavepacket
    velocity = 0.5  # Average velocity of the packet
    k0 = mass*velocity/h_bar       # Average wavenumber
    sigma0 = L/10   # Standard deviation of the wavefunction
    Norm_psi = 1/(np.sqrt(sigma0*np.sqrt(np.pi)))  # Normalization
    psi = Norm_psi * np.exp(i_imag*k0*x) *np.exp(-(x-x0)**2/(2*sigma0**2))
    
    E = (h_bar*k0)**2/(2*mass)
    print(f"Particle Energy is: {E}")
    
    #* Set up the Hamiltonian operator matrix
    ham = np.zeros((N,N))  # Set all elements to zero
    coeff = -h_bar**2/(2*mass*h**2)
    for i in range(1,(N-1)):
      ham[i,i-1] = coeff
      ham[i,i] = -2*coeff  # Set interior rows
      ham[i,i+1] = coeff
    
    # First and last rows for periodic boundary conditions
    ham[0,N-1] = coeff
    ham[0,0] = -2*coeff
    ham[0,1] = coeff
    ham[N-1,N-2] = coeff
    ham[N-1,N-1] = -2*coeff
    ham[N-1,0] = coeff
    
    delta = 1
    if delta:
        #U = float(input("Enter delta function amplitude:")) Take U from function args
        # Now we add the delta function right in the middle of our interval
        """
        # If odd number of points we put the delta func only on the middle point
        if N%2 == 1:
            ham[int(N/2),int(N/2)] += U*(1/h) # U times magnitude of delta function
        elif N%2 == 0: # If even number of grid points, we put the delta func on the middle 2 points
            ham[int(N/2)+1,int(N/2)+1] += U*(1/(2*h))
            ham[int(N/2)+1,int(N/2)+1] += U*(1/(2*h))     
            
            """
                # If odd number of points we put the delta func only on the middle point
        if N%2 == 1:
            ham[:,int(N/2)] += U*(1/h) # U times magnitude of delta function
        elif N%2 == 0: # If even number of grid points, we put the delta func on the middle 2 points
            ham[:,int(N/2)+1] += U*(1/(2*h))
            ham[:,int(N/2)+1] += U*(1/(2*h))  
    #* Compute the Crank-Nicolson matrix
    dCN = np.dot( np.linalg.inv(np.identity(N) + 0.5*i_imag*tau/h_bar*ham),(np.identity(N) - .5*i_imag*tau/h_bar*ham))
    
    
    
    #* Plot the initial wavefunction
    
    plt.figure(1); plt.clf() 
    plt.plot(x,np.real(psi),'-',x,np.imag(psi),'--')
    plt.title('Initial wave function')
    plt.xlabel('x')
    plt.ylabel(r'$\psi(x)$')
    plt.legend(['Real  ','Imag  '])
    
    # pause(1)
    
    plt.figure(2)
    # plt.ion()
    #* Initialize loop and plot variables 
    max_iter = int(L/(velocity*tau) ) +1    # Particle should circle system
    #max_iter = int(max_iter/10) # Plot 1/10th of a rotation
    print(f'Max_iter: {max_iter}')
    plot_iter = int(max_iter/5)         # Produce 10 curves modded to 5
    # For some reason, I have to convert to a real number, even though a complex number
    # times its conjugate should be real, numpy still thinks it is complex and gets upset
    p_plot = np.real(psi*np.conj(psi) )    # Record initial condition
    iplot = 0
    tplot = np.array([0])
    axisV = np.array([-L/2., L/2., 0., np.max(p_plot)])# Fix axis min and max
    
    #* Loop over desired number of steps (wave circles system once)
    for iter in range(1,max_iter):
    	
      #* Compute new wave function using the Crank-Nicolson scheme
      psi = np.dot(dCN,psi)
      
      #* Periodically record values for plotting
      if(iter%plot_iter == 0 ):
        iplot = iplot+1
        p_plot = np.vstack((p_plot,np.real(psi*np.conj(psi)) ))
        tplot = np.append(tplot,tau*iter)
    
    #* Plot probability versus position at various times
    
    plt.figure(3); plt.clf() 
    pFinal = np.real(psi*np.conj(psi))
    for i in range(iplot):
        plt.plot(x,p_plot[i,:],label='t='+str(tplot[i]))
    plt.plot(x,pFinal)
    plt.xlabel('x')
    plt.ylabel('P(x,t)')
    plt.legend()
    plt.title('Probability density at various times')
    plt.show()
   
    return