# -*- coding: utf-8 -*-
"""
Phys 416
HW 1
Alexander Swart
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def graph():
    
    #Define x-axes
    x = np.linspace(0,20,1000)
    sub_x = np.linspace(0,20,6)
    theta = np.linspace(0,2*np.pi,1000)
    
    
    # Define Subplots
    fig,ax = plt.subplots(2,2)
    ax[1,1] = plt.subplot(224,polar=True) # One must be polar
    
    # Define Equations
    y1 = np.exp(-x/4)*np.sin(x)
    sub_y1 = np.exp(-sub_x/4)*np.sin(sub_x)
    y2 = np.exp(-x)*(np.sin(x)**2)
    y3 = np.exp(-x/4)
    
    r1 = np.cos(theta*3)**2
    
    # Plot x,y coordinates
    ax[0,0].plot(x,y1)
    
    ax[0,1].plot(x,y1,color='C0')
    ax[0,1].plot(x,y3,c='red',linestyle='--')
    #ax[0,1].scatter(sub_x,sub_y1)
    ax[0,1].plot(sub_x,sub_y1,color='green',marker='o',linestyle='')

    
    ax[1,0].scatter(x,y2,s=8,c='k')
    ax[1,1].plot(theta,r1)
    
    # Formatting
    fig.tight_layout(pad=2,h_pad=2,w_pad=.5)
    

    
    # Set titles
    ax[0,0].set_title(r'$y = e^{-\frac{x}{4}}\sin{(x)}$')
    ax[0,1].set_title(r'$y = e^{-\frac{x}{4}}\sin{(x)}$')
    ax[1,0].set_title(r'$y = e^{-x}\sin^{2}{(x)}$')
    
    # Set x-axis labels
    ax[0,0].set_xlabel('x')
    ax[0,1].set_xlabel('x')
        
    # Set axis limits
    ax[0,0].set_ylim(-.5,1)
    ax[0,1].set_ylim(-.5,1)    
    ax[1,1].set_ylim(0,1.05)
    
    # Set log axis scale
    ax[1,0].set_yscale('log')
    ax[1,0].set_ylim(ymin=10**-10,ymax=10**0)

    # Set interior plot margins
    ax[0,0].margins(0)
    ax[0,1].margins(0)
    ax[1,0].margins(0)

    # Set axis ticks
    cartesian_xticks = [0,5,10,15,20]
    ax[0,0].set_xticks(cartesian_xticks)
    ax[0,1].set_xticks(cartesian_xticks)
    ax[1,0].set_xticks(cartesian_xticks)
    
    ax[0,0].set_yticks([-.5,0,.5,1])
    ax[0,1].set_yticks([-.5,0,.5,1])
    ax[1,0].set_yticks([10**-10,10**-5,10**0])
    ax[1,0].get_yaxis().get_major_formatter().labelOnlyBase = False # Setting log axis ticks requires this
    ax[1,1].set_yticks([0,0.5,1])
    
    # Set grids on axis ticks
    ax[0,0].grid(True,color='b',linestyle='--')
    ax[0,1].grid(True,color='b',linestyle='--')
    ax[1,0].grid(True,color='b',linestyle='--')
    ax[1,1].grid(True)
    

    return

if __name__ == '__main__':
    plt.close()
    graph()
    