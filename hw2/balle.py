# python 3 version
import matplotlib.pyplot as plt
import numpy as np
# balle - Program to compute the trajectory of a baseball
#         using the Euler method.


def balle(tau = .1, get_input = False, calc_error = False, plot_energy = True, midpoint = False):
    
    
    # Get input values from input prompts
    if get_input:
    #* Set initial position and velocity of the baseball  
        y1 = float(input("Enter initial height (meters): "))
        speed = float(input("Enter initial speed (m/s): "))
        theta = float(input("Enter initial angle (degrees): "))
        airFlag = bool(input("Air resistance? (Yes:1, No:0):"))
        tau = float(input("Enter timestep, tau (sec): "));  # (sec)
    else:
        # Set default initial conditions for experimenting with tau
        y1 = 0.0
        speed = 5.0
        theta = 45.0
        airFlag = False
        
        
    r1 = np.array([0.0, y1]);     # Initial vector position    
    v1 = np.array([speed*np.cos(theta*np.pi/180), speed*np.sin(theta*np.pi/180)])  # Initial velocity
    r = np.copy(r1)
    v = np.copy(v1)  # Set initial position and velocity, best to copy to avoid overwrites
    
    #* Set physical parameters (mass, Cd, etc.)
    Cd   = 0.35;    # Drag coefficient (dimensionless)
    area = 4.3e-3;  # Cross-sectional area of projectile (m^2)
    grav = 9.81;    # Gravitational acceleration (m/s^2)
    mass = 0.145;   # Mass of projectile (kg)
    
    
    if airFlag:
        rho = 0;      # No air resistance
    else: 
        rho = 1.2;    # Density of air (kg/m^3)
    
    air_const = -0.5*Cd*rho*area/mass;  # Air resistance constant
    
    #* Loop until ball hits ground or max steps completed
    maxstep = 10000;   # Maximum number of steps
    for istep in range(0,maxstep):
        #* Record position (computed and theoretical) for plotting
        t = (istep)*tau;     # Current time
        if(istep ==0):
            xplot = np.array(r[0]);   # Record trajectory for plot
            yplot = np.array(r[1]);
            xNoAir = np.array(r[0]);
            yNoAir = np.array(r[1]);
            time = np.array(t)
            velocity = np.array(v)
        else:
            xplot = np.append(xplot,r[0]);   # Record trajectory for plot
            yplot = np.append(yplot,r[1]);
            xNoAir = np.append(xNoAir,r1[0] + v1[0]*t);   # Record trajectory for plot
            yNoAir = np.append(yNoAir,r1[1] + v1[1]*t - 0.5*grav*t**2);   
             
        #* Calculate the acceleration of the ball 
        accel = air_const*np.linalg.norm(v)*v;   # Air resistance
        accel[1] = accel[1]-grav;      # Gravity
    
        #* Calculate the new position and velocity using Euler method
       
        #v_new = v + tau*accel # Midpoint method
        #r = r + (tau/2)*(v+v_new)  #Midpoint method
        r = r + (tau)*(v)                 # Euler step
        v = v + tau*accel
        #v = v_new #Midpoitn method
        time = np.append(time,t)
        velocity = np.append(velocity,v)
        #* If ball reaches ground (y<0), break out of the loop
        if( r[1] < 0 ):  
            xplot = np.append(xplot,r[0]);   # Record trajectory for plot
            yplot = np.append(yplot,r[1]);
            time = np.append(time,t*tau)
            break;                  # Break out of the for loop
    
    # Once the ball reaches the ground, interpolate the last 3 points to find accurate endpoints
    x_end = np.interp(0,yplot[-3:],xplot[-3:]) # Note use interpf
    t_end = np.interp(0,yplot[-3:],time[-3:])
    print(x_end,t_end)
    
    
    # Print maximum range and time of flight
    print('Maximum range is ',r[0],' meters');
    print('Time of flight is ',istep*tau,' seconds');
    
    # Graph the trajectory of the baseball
    plt.figure(0);   plt.clf(); # Clear figure window and bring it forward
    # Mark the location of the ground by a straight line
    xground = np.array([0, np.max(xNoAir)]);  yground = np.array([0, 0]);
    # Plot the computed trajectory and parabolic, no-air curve
    plt.plot(xplot,yplot,'.')
    plt.plot(xNoAir,yNoAir,'--');
    plt.plot(xground,yground,'-');
    plt.legend(['Euler method','Theory (No air)']);
    plt.xlabel('Range (m)');  plt.ylabel('Height (m)');
    plt.title('Projectile motion');
    #axis equal; shg; # reset the aspect ratio, bring the plot to the front
    plt.grid(True)
    plt.show()
    
    if plot_energy:
        plt.figure(1)
        
    
    
    
    return velocity
    


if __name__ == '__main__':
    v = balle()