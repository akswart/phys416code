"""
FEniCS tutorial demo program: Incompressible Navier-Stokes equations
for channel flow (Poisseuille) on the unit square using the
Incremental Pressure Correction Scheme (IPCS).

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0

See 1st part of: https://fenicsproject.org/pub/tutorial/html/._ftut1009.html
"""
import timeit
import multiprocessing
import matplotlib.pyplot as plt
from fenics import *
import numpy as np
from tqdm import tqdm # Way lower overhead than ProgressBar 60ns per iter vs 800ns see github



def navier(grid_size):
    set_log_level(LogLevel.ERROR)

    T = 10.0           # final time
    num_steps = 500    # number of time steps
    dt = T / num_steps # time step size
    mu = 1             # kinematic viscosity
    rho = 1            # density
    
    # Create mesh and define function spaces
    mesh = UnitSquareMesh(grid_size, grid_size)
    V = VectorFunctionSpace(mesh, 'P', 2)
    Q = FunctionSpace(mesh, 'P', 1)
    
    # Define boundaries
    inflow  = 'near(x[0], 0)'
    outflow = 'near(x[0], 1)'
    walls   = 'near(x[1], 0) || near(x[1], 1)'
    
    # Define boundary conditions
    bcu_noslip  = DirichletBC(V, Constant((0, 0)), walls)
    bcp_inflow  = DirichletBC(Q, Constant(8), inflow)
    bcp_outflow = DirichletBC(Q, Constant(0), outflow)
    bcu = [bcu_noslip]
    bcp = [bcp_inflow, bcp_outflow]
    
    # Define trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)
    
    # Define functions for solutions at previous and current time steps
    u_n = Function(V)
    u_  = Function(V)
    p_n = Function(Q)
    p_  = Function(Q)
    
    # Define expressions used in variational forms
    U   = 0.5*(u_n + u)
    n   = FacetNormal(mesh)
    f   = Constant((0, 0))
    k   = Constant(dt)
    mu  = Constant(mu)
    rho = Constant(rho)
    
    # Define strain-rate tensor
    def epsilon(u):
        return sym(nabla_grad(u))
    
    # Define stress tensor
    def sigma(u, p):
        return 2*mu*epsilon(u) - p*Identity(len(u))
    
    # Define variational problem for step 1
    F1 = rho*dot((u - u_n) / k, v)*dx + \
         rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
       + inner(sigma(U, p_n), epsilon(v))*dx \
       + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
       - dot(f, v)*dx
    a1 = lhs(F1)
    L1 = rhs(F1)
    
    # Define variational problem for step 2
    a2 = dot(nabla_grad(p), nabla_grad(q))*dx
    L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx
    
    # Define variational problem for step 3
    a3 = dot(u, v)*dx
    L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx
    
    # Assemble matrices
    A1 = assemble(a1)
    A2 = assemble(a2)
    A3 = assemble(a3)
    
    # Apply boundary conditions to matrices
    [bc.apply(A1) for bc in bcu]
    [bc.apply(A2) for bc in bcp]
    
    # Time-stepping
    t = 0
    time_array = np.linspace(t,T,num_steps)
    error_array = np.zeros(num_steps)
    for n in range(num_steps):
    
        # Update current time
        t += dt
    
        # Step 1: Tentative velocity step
        b1 = assemble(L1)
        [bc.apply(b1) for bc in bcu]
        solve(A1, u_.vector(), b1)
    
        # Step 2: Pressure correction step
        b2 = assemble(L2)
        [bc.apply(b2) for bc in bcp]
        solve(A2, p_.vector(), b2)
    
        # Step 3: Velocity correction step
        b3 = assemble(L3)
        solve(A3, u_.vector(), b3)
    
        # Plot solution
        plot(u_)
    
        # Compute error
        u_e = Expression(('4*x[1]*(1.0 - x[1])', '0'), degree=2)
        u_e = interpolate(u_e, V)
        error = np.abs(np.array(u_e.vector()) - np.array(u_.vector())).max()
        error_array[n] = error
        #print('t = %.2f: error = %.3g' % (t, error))
        #print('max u:', np.array(u_.vector()).max())
    
        # Update previous solution
        u_n.assign(u_)
        p_n.assign(p_)
    
    return error_array,time_array

def navier_benchmark(grid_size):
    # Run once to get error_array, time_array
    error_array,time_array = navier(grid_size)
    
    # Run with timeit to get average time
    mysetup = f"""grid_size = {grid_size}
from ft07_navier_stokes_channel import navier                
    """
    mystmt = "navier(grid_size)"
    t = timeit.timeit(setup = mysetup,stmt = mystmt,number=10)
    return (t,error_array,time_array)


if __name__ == "__main__":

    grid_size = np.arange(2,17,2)

    p = multiprocessing.Pool()

    result = p.map(navier_benchmark,grid_size)
    mapped_results = {grid_size[i]:result[i] for i in range(len(grid_size))}
    p.close()
    p.join()

    times =  np.array([i[0] for i in result])
    errors = np.array([i[1] for i in result])
    time_arrays = np.array([i[2] for i in result])
    np.savez("benchmarkout.npz",grid_size,times,errors,time_arrays)
    plt.figure()
    plt.scatter(grid_size,times)
    plt.show()
    
    