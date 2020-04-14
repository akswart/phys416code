# non Neumann demo - ftcs

import numpy as np
import matplotlib.pyplot as plt

angle = np.arange(0,2*np.pi,np.pi/50.)

d=np.array([0.125, 0.25, 0.5, 0.6 ]) # d = Tau * kappa/h^2
r1=np.ones(len(angle))

# amplification factor as a function of d and angle
a=np.zeros([len(angle),len(d)])
aa=np.copy(a)
x=np.zeros([len(angle),len(d)])
y=np.copy(x)
for j in range(0,len(d)):
    for i in range(0,len(angle)):
    
        #a[i,j]=1-2*d[j]*(1-np.cos(angle[i])) Orig a
        # Modified a
        a[i,j] = 2*d[j]*(np.cos(angle[i])-1) - np.sqrt( (2*d[j])**2*(np.cos(angle[i]) - 1)**2 + 1)
        aa[i,j]=np.abs(a[i,j])
#  for plotting
        x[i,j] = a[i,j]*np.cos(angle[i])
        y[i,j] = a[i,j]*np.sin(angle[i])



print(str(d[0:]))
plt.figure(1)
plt.clf()
plt.plot(x[:,0],y[:,0],x[:,1],y[:,1],x[:,2],y[:,2],x[:,3],y[:,3])
leg=['d='+str(d[0]),' d='+str(d[1]),' d='+str(d[2]),' d='+str(d[3])]
plt.grid(True)
plt.legend(leg)
plt.axis('equal')
plt.plot(np.cos(angle),np.sin(angle),'r--',linewidth=2)
plt.grid(True)
plt.title('A')

plt.figure(2)
plt.clf()
plt.plot(angle,a[:,0],angle,a[:,1],angle,a[:,2],angle,a[:,3])
plt.grid(True)
plt.ylabel(' A')
plt.xlabel('angle kh (radians)')
plt.legend(leg)

plt.show()