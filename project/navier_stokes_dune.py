"""
Incompressible Navier-Stokes equations
for flow over a sand dune using the Incremental Pressure Correction
Scheme (IPCS).

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0

See https://fenicsproject.org/pub/tutorial/html/._ftut1009.html for theory and 
unmodified code.

Also, for some wierd reason, its wasnt writing correctly to the hdf5 and xdmf 
outputs when I ran it in spyder, but ran perfectly from command line. No idea why.

"""

import fenics as fs
from fenics import dot,dx,nabla_grad,ds,inner,div # Import math operations by name to clean up code
import mshr
from dune_mesh import generate_dune_mesh
#from mshr import *
import numpy as np
import time

from tqdm import tqdm # Way lower overhead than ProgressBar 60ns per iter vs 800ns see github

import os
t0 = time.time()


# Clean up output directory before running

if os.name == 'nt': #For windows systtems
    dirs = os.listdir()
    if 'navier_stokes_dune' in dirs:
        try: # Do a try incase zipping fails
            #We will append this time to the filename to make it unique 
            sys_time = time.time() # Num of seconds since the epoch

            from zipfile import ZipFile 
  
            def get_all_file_paths(directory): 
              
                # initializing empty file paths list 
                file_paths = [] 
              
                # crawling through directory and subdirectories 
                for root, directories, files in os.walk(directory): 
                    for filename in files: 
                        # join the two strings in order to form the full filepath. 
                        filepath = os.path.join(root, filename) 
                        file_paths.append(filepath) 
              
                # returning all file paths 
                return file_paths  

                # path to folder which needs to be zipped 
            directory = './navier_stokes_dune'
          
            # calling function to get all file paths in the directory 
            file_paths = get_all_file_paths(directory) 
          
            # printing the list of all files to be zipped 
            print('Following files will be zipped:') 
            for file_name in file_paths: 
                print(file_name) 
          
            # writing files to a zipfile
            print("Zipping old output files")
            with ZipFile(f'navier_stokes_dune{sys_time}.zip','w') as zip: 
                # writing each file one by one 
                for file in tqdm(file_paths): 
                    zip.write(file) 
            print('All files zipped successfully!')
            os.system("rmdir /Q /S navier_stokes_dune") # Remove old folder
        except:
            raise UserWarning("compressed backup not made, output dir will be deleted")
            os.system('rm -rf navier_stokes_dune')

elif os.name == 'posix': # For linux and mac
    dirs = os.listdir()
    if "archive" not in dirs: # If we dont already have an archive, make one
        os.system("mkdir archive")
        print("Made archive folder to store old results")
    else:
        print("Putting old output into folder: archive")
    if 'navier_stokes_dune' in dirs:
        try: # Do a try, except in case tar is not installed to zip files
            #We will append this time to the filename to make it unique 
            sys_time = time.time() # Num of seconds since the epoch
            print("Compressing old files")
            os.system(f"tar -czvf ./archive/navier_stokes_dune{sys_time}.tar.gz navier_stokes_dune")
            print("Deleting old files")
            os.system('rm -rf navier_stokes_dune')
        except: # If tar is not installed attempt to do it with zipfile but its not as fast as tar
            #raise UserWarning("tar not found, compressed backup not made, output dir will be deleted")
            #os.system('rm -rf navier_stokes_cylinder')
            #"""
            raise UserWarning("tar not found")
            try:
                #We will append this time to the filename to make it unique 
                sys_time = time.time() # Num of seconds since the epoch
    
                from zipfile import ZipFile 
      
                def get_all_file_paths(directory): 
                  
                    # initializing empty file paths list 
                    file_paths = [] 
                  
                    # crawling through directory and subdirectories 
                    for root, directories, files in os.walk(directory): 
                        for filename in files: 
                            # join the two strings in order to form the full filepath. 
                            filepath = os.path.join(root, filename) 
                            file_paths.append(filepath) 
                  
                    # returning all file paths 
                    return file_paths  
    
                # path to folder which needs to be zipped 
                directory = './navier_stokes_dune'
              
                # calling function to get all file paths in the directory 
                file_paths = get_all_file_paths(directory) 
              
                # printing the list of all files to be zipped 
                print('Following files will be zipped:') 
                for file_name in file_paths: 
                    print(file_name) 
              
                # writing files to a zipfile
                print("Zipping old output files")
                with ZipFile(f'navier_stokes_dune{sys_time}.zip','w') as zip: 
                    # writing each file one by one 
                    for file in tqdm(file_paths): 
                        zip.write(file) 
                print('All files zipped successfully!')                
                os.system('rm -rf navier_stokes_dune')    
            except:
                raise UserWarning("tar not found, compressed backup not made, output dir will be deleted")
                os.system('rm -rf navier_stokes_dune')
   
else: # System not recognized
    raise UserWarning("system not recognized, compressed backup not made, output dir must be deleted manually")


print("Starting to solve incompressible Navier–Stokes equations")
T = 5.0            # final time
num_steps = 10000   # number of time steps
dt = T / num_steps # time step size
print(f"Timestep: {dt}")
mu = 0.001         # dynamic viscosity
rho = 1            # density

print("Creating mesh")
# Create mesh, from dune_mesh
large_scale = True
mesh = generate_dune_mesh(large_scale,64)  # True is large scale view, using 16 grid points initaly for fast testing runs
#mesh = mshr.generate_mesh(mshr.Rectangle(fs.Point(0, 0), fs.Point(2.2, .41)),16)

# Define function spaces
V = fs.VectorFunctionSpace(mesh, 'P', 2)
Q = fs.FunctionSpace(mesh, 'P', 1)


print("Defining boundary conds")
# Define boundaries
if large_scale: # Different boundaries for zoomed out view x=(-4,20) y=(0,10)
    inflow   = 'near(x[0], 0)'
    outflow  = 'near(x[0], 2.2)'
    walls    = 'near(x[1], .41)'
    # The dune goes from x(-2,10) y=(0,1)
    # So, if we capture all points in this region and part of a boundary, 
    # we can set the boundary conditions for the dune
    dune = 'near(x[1], 0) || on_boundary && x[0]>.02 && x[0]<2 && x[1]>=0 && x[1]<.3'
    
    
#inflow   = 'near(x[0], 0)' # Inflow and outflow are unchanged
#outflow  = 'near(x[0], 2.2)'
#walls    = 'near(x[1], 0) || near(x[1], 0.41)'
#cylinder = 'on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3'


# Define inflow profile
flow_type = 2

if flow_type == 1:
    inflow_profile = ('4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0') # profiled wind as in Fenics example
elif flow_type == 2:
    inflow_profile = ('1', '0') # uniform wind


# Define boundary conditions
bcu_walls = fs.DirichletBC(V, fs.Constant((0, 0)), walls)
bcu_dune = fs.DirichletBC(V, fs.Constant((0, 0)), dune)

bcu_inflow = fs.DirichletBC(V, fs.Expression(inflow_profile, degree=2), inflow)

#bcp_inflow = fs.DirichletBC(Q, fs.Constant(.5), inflow)
bcp_outflow = fs.DirichletBC(Q, fs.Constant(0), outflow)
bcu = [bcu_inflow,bcu_walls, bcu_dune]
#bcu = [bcu_inflow,bcu_walls]
bcp = [ bcp_outflow]

# Define trial and test functions
u = fs.TrialFunction(V)
v = fs.TestFunction(V)
p = fs.TrialFunction(Q)
q = fs.TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = fs.Function(V)
u_  = fs.Function(V)
p_n = fs.Function(Q)
p_  = fs.Function(Q)

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n  = fs.FacetNormal(mesh)
f  = fs.Constant((0, 0))
k  = fs.Constant(dt)
mu = fs.Constant(mu)
rho = fs.Constant(rho)

# Define symmetric gradient
def epsilon(u):
    return fs.sym(fs.nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*fs.Identity(len(u))


# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx \
   + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
   - dot(f, v)*dx
a1 = fs.lhs(F1)
L1 = fs.rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = fs.assemble(a1)
A2 = fs.assemble(a2)
A3 = fs.assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Create XDMF files for visualization output
xdmffile_u = fs.XDMFFile('navier_stokes_dune/velocity.xdmf')
xdmffile_p = fs.XDMFFile('navier_stokes_dune/pressure.xdmf')

hdf = fs.HDF5File(mesh.mpi_comm(), "navier_stokes_dune/file.h5", "w")
hdf.write(mesh, "/mesh")


# Create time series (for use in reaction_system.py)
timeseries_u = fs.TimeSeries('navier_stokes_dune/velocity_series')
timeseries_p = fs.TimeSeries('navier_stokes_dune/pressure_series')

# Save mesh to file (for use in reaction_system.py)
#File('navier_stokes_cylinder/cylinder.xml.gz') << mesh

# Create progress bar
progress = fs.Progress('Time-stepping',num_steps)
fs.set_log_level(fs.LogLevel.PROGRESS)


# Time-stepping
print("Time Stepping")
t = 0
for n in tqdm(range(num_steps)):

    fs.set_log_level(fs.LogLevel.ERROR) # Only log explody stuff
    # Update current time
    t += dt

    # Step 1: Tentative velocity step
    b1 = fs.assemble(L1)
    [bc.apply(b1) for bc in bcu]
    fs.solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')

    # Step 2: Pressure correction step
    b2 = fs.assemble(L2)
    [bc.apply(b2) for bc in bcp]
    fs.solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')

    # Step 3: Velocity correction step
    b3 = fs.assemble(L3)
    fs.solve(A3, u_.vector(), b3, 'cg', 'sor')

    # Plot solution
    #fs.plot(u_, title='Velocity')
    #fs.plot(p_, title='Pressure')

    # Save solution to file (XDMF/HDF5)
    xdmffile_u.write(u_, t)
    xdmffile_p.write(p_, t)

    # Save nodal values to file
    timeseries_u.store(u_.vector(), t)
    timeseries_p.store(p_.vector(), t)

    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)

    # Update progress bar
    #set_log_level(LogLevel.PROGRESS)
    progress += 1
    #print('u max:', u_.vector().array().max())
    
tf = time.time()
print(tf-t0)
