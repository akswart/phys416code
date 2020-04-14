# advect - Program to solve the advection equation 
# using the various hyperbolic PDE schemes - high resolution version
# for questions 1 and 2
# clear all  help advecth  # Clear memory and print header
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def hires2(y,u,h,tau,limit):
# hi resolution function that uses a limiter

    N=len(y)
    n2=len(u)

    if N != n2:
        print(' lengths of y and u do not match')
        return y
#        yout(1:N) = y(1:N) - tau/h*u(1:N).*(y(1:N)-y(im));

    uplus =np.maximum(u,0)
    uminus=np.minimum(u,0)

    delta = np.zeros(N)
    yout = np.zeros(N)
    delta = limiter2(y,u,limit)

    for i in range(0,N):
        flux=uminus[i]*y[i]+uplus[i]*y[im[i]]+0.5*abs(u[i])*(1-abs(u[i]*tau/h))*delta[i]
        fluxp=uminus[i]*y[ip[i]]+uplus[i]*y[i]+0.5*abs(u[i])*(1-abs(u[i]*tau/h))*delta[ip[i]]
        yout[i] = y[i] -tau/h*(fluxp-flux)

    return yout

def limiter2(y,u,limit):
    N=len(y)
    n2=len(u)
    if N != n2:
        print(' lengths of y and u do not match')
        return y

    deltay=y[0:N]-y[im]
    I=np.copy(ip)
    ipositive = np.where(u>0)
    I[ipositive]=im[ipositive]

    theta = np.zeros(N)
    for i in range(0,N):
        theta[i]=0.0
        if deltay[i] != 0:
            theta[i]=deltay[I[i]]/deltay[i]

#limit = 'minmod'
    phi = np.zeros(N)
    if(limit=='upwind'):
        phi=np.zeros(N) # upwind
    elif(limit=='lw'):
        phi = np.ones(N) # lax-wendroff
    elif(limit=='bm'):
        phi=np.copy(theta) # beam warming
    elif(limit=='minmod'):
        phi[0:N]=minmod(np.ones(N),theta[0:N]) # minmod method
    elif(limit =='mc'):
        for i in range(0,N):
          phi[i] = np.max([0,np.min([(1+theta[i])/2,2,2*theta[i]])]) # MC limiter
    elif(limit =='vanleer'):
         phi[0:N] = (theta[0:N] + abs(theta[0:N]))/(1+abs(theta[0:N])) # van leer
    elif(limit =='superbee'):
        phi[:] = 0.0 # superbee limiter
    else:
        print('Unknown limiter method.')

    delta=phi*deltay
    return delta
def minmod(a,b):
# minmod function - array smart

    if len(a) != len(b):
        print(' minmod, sizes do not match')

    phi = ((a*b)>0)*((a*(abs(a)<=abs(b)))+(b*(abs(a)>abs(b))))
    return phi
#* Select numerical parameters (time step, grid spacing, etc.).

method = int(input('Choose a numerical method: 1- FTCS, 2-Lax, 3-Lax-Wendroff, 4-Upwind, 5-High res: '))

if method ==5: # select limiter
    choice = int(input('For the hires method, choose a limiter: 1-upwind,2-Lax-Wendroff, \
    3-beam warming, 4-minmod, 5-MC, 6-van leer, 7-superbee: '))
    # convert to text-based limit to pass to the limiter
    if (choice==1):
        limit='uw'
    elif (choice==2):
        limit='lw'
    elif (choice==3):
        limit='bm'
    elif (choice==4):
        limit='minmod'
    elif (choice==5):
        limit='mc'
    elif (choice==6):
        limit='vanleer'
    elif (choice==7):
        limit='superbee'

N = int(input('Enter number of grid points: '))
L = 1.     # System size
h = L/N    # Grid spacing
c = 1      # Wave speed
print('Time for wave to move one grid spacing is ',(h/c))
tau = float(input('Enter time step: '))
coeff = -c*tau/(2.*h)  # Coefficient used by all schemes
coefflw = 2*coeff**2    # Coefficient used by L-W scheme
print('Wave circles system in %d steps'%(L/(c*tau)))
nStep = int(input('Enter number of steps: '))

#* Set initial and boundary conditions.
sigma = 0.1              # Width of the Gaussian pulse
k_wave = np.pi/sigma        # Wave number of the cosine
x = (np.arange(0,N)+1/2)*h - L/2  # Coordinates of grid points

ic=int(input('Input initial condition:, 1-gaussian pulse, 2-square wave, 3-both: '))
#* Set initial and boundary conditions.
if(ic ==1):
    sigma = 0.1              # Width of the Gaussian pulse
    k_wave = np.pi/sigma        # Wave number of the cosine
    # Initial condition is a Gaussian-cosine pulse
    a = np.cos(k_wave*x) * np.exp(-x**2/(2*sigma**2))
elif(ic==2):
    a=np.zeros(N)
    for i in range(int(N/4),int(N/2)):
        a[i]= 1.
else:
    sigma = 0.025;              # Width of the Gaussian pulse
    k_wave = np.pi/sigma;        # Wave number of the cosine
    # Initial condition is a Gaussian-cosine pulse
    a =  np.exp(-(x-L/4)**2/(2*sigma**2))
    for i in range(int(N/4),int(N/2)):
        a[i] = 1.0


# Use periodic boundary conditions
# Use periodic boundary conditions
ip = np.arange(0,N)+1
ip[N-1] = 0 # ip = i+1 with periodic b.c.
im = np.arange(0,N)-1
im[0] = N-1 # im = i-1 with periodic b.c.

#* Initialize plotting variables.
iplot = 0          # Plot counter
aplot = np.copy(a)  # Record the initial state
tplot = np.array([0])       # Record the initial time (t=0)
nplots = 50        # Desired number of plots
plotStep = max(1, np.floor(nStep/nplots)) # Number of steps between plots

#* Loop over desired number of steps.
# plt.ion() # this messes things up
for iStep in range(nStep+1):  ## MAIN LOOP ##
    
    #* Compute new values of wave amplitude using FTCS, 
    #  Lax or Lax-Wendroff method.
    if( method == 1 ):     ### FTCS method ###
        a[0:N]= a[0:N] + coeff*(a[ip]-a[im])
    elif( method == 2 ):  ### Lax method ###
        a[0:N] = 0.5*(a[ip]+a[im]) + coeff*(a[ip]-a[im])
    elif( method==3):                 ### Lax-Wendroff method ###
        a[0:N] = a[0:N] + coeff*(a[ip]-a[im]) + coefflw*(a[ip]+a[im]-2*a)
    elif( method==4):                ### Upwind method ###
        a[0:N] = a[0:N] + 2*coeff*(a[0:N]-a[im])
    elif( method==5):   ### Hi res method ###
        u=c*np.ones(len(a))
        a = hires2(a,u,h,tau,limit)
    
    #* Periodically record a(t) for plotting.

    if( (iStep%plotStep) < 1):  # Every plot_iter steps record
        iplot = iplot+1
        aplot = np.vstack((aplot,a))       # Record a(i) for plotting
        tplot = np.append(tplot,tau*iStep)
        print('%d out of %d steps completed'%(iStep,nStep))
        


#* Plot the initial and final states.
# need plt.ion() for plot windows to update
animate=1
if(animate==1):
    # plots in a movie fashion - comment out if you want the program to be faster
    # plt.ion()    
    for i in range(iplot+1):
        plt.figure(1)   # Clear figure 1 window and bring forward
        plt.clf()
        plt.plot(x,aplot[0,:],'-',label='Initial')
        plt.plot(x,aplot[i,:],'--',label='current')
        plt.legend(['Initial  ','Final'])
        plt.xlabel('x')
        plt.ylabel('a(x,t)')
        plt.grid(True)
        plt.axis([-0.5, 0.5, -0.5, 1.2])
        plt.title(' time ='+str(tplot[iplot]))
        if(method == 1):
            plt.title('FTCS method, time ='+str(tplot[i]))
        elif(method == 2):
            plt.title('Lax method, time =' +str(tplot[i]))
        elif(method == 3):
            plt.title('Lax-Wendroff method, time =' +str(tplot[i]))
        elif(method == 4):
            plt.title('Upwind method, time='+str(tplot[i]))
        elif(method==5):
            plt.title('High resolution method, time ='+str(tplot[i]))
    
        plt.legend()
        plt.draw()
        # if iStep == nStep-1:
        #     temp=input('Hit any key to stop')
        plt.pause(tau)
        # end plots in a movie fashion
        
        plt.show()
    # plt.ioff()

plt.figure(2)
plt.clf()   # Clear figure 2 window and bring forward
plt.plot(x,aplot[0,:],'-',x,a,'--')
plt.legend(['Initial  ','Final'])
plt.xlabel('x')
plt.ylabel('a(x,t)')
plt.axis([-0.5, 0.5, -0.5, 1.2])
plt.grid(True)
if (method == 1):
    plt.title('FTCS method')
elif(method == 2):
    plt.title('Lax method')
elif(method == 3):
    plt.title('Lax-Wendroff method')
elif(method == 4):
    plt.title('Upwind method')
elif(method==5):
    plt.title('High resolution using the '+ limit + ' method')
plt.show()

# #* Plot the wave amplitude versus position and time
tt,xx = np.meshgrid(x,tplot)
fig = plt.figure(3);plt.clf()   # Clear figure 3 window and bring forward
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xx, tt, aplot, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
ax.set_xlabel('Time')
ax.set_ylabel('Position')
ax.set_zlabel('Amplitude)')

plt.show()
