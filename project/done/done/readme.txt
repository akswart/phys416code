readme

FOR FURTHER EXPLANATION I HIGHLY ENCOURAGE YOU TO READ THE JUPYTER NOTEBOOK, Dune_sim.ipynb

The following will walk you through on setting up the right environment and
running the example code provided.

install anaconda
Create the environment from the environment.yml file with:
>    conda env create -f environment.yml
Activate the new environment
>    conda activate dunephys
Run the example file:
>    python navier_stokes_dune_small_run.py
This will run an example of the navier_stokes_done.py code
with variables appropriate for local execution
Note: In order to achive the resolution as in the output animations, 
navier_stokes_dune.py must be run on a larger grid. This can take hours.

It should take about 1 min on a reasonably modern laptop, and 
will use parameters T=2.0, num_steps = 2000, and a mesh_size of 32

This will produce a folder called navier_stokes_dune in which 
the results will be stored. In order to visualize these, you need something capable of viewing hdf5 files, I used Paraview (another option is VisIT)

In paraview, in order to view the velocity, for example, you will open up the 
velocity.xdmf file with a compatible reader, then hit apply, then go down to the Coloring section in the Properties tab, select the box next to "Magnitude that should say "Solid Color" then select f_31. Finally hit run and you will be able to see the valocity of the fluid change in real time.
