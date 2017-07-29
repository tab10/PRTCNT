[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://github.com/tab10/conduction/LICENSE)

# PRTCNT
## (Python)(Randomwalksimulationforobtaining)(Thermalconductiviesof)(Carbon)(Nano)(Tubecomposites)

#### Python/MPI-based simulation package to measure the thermal transport properties of a carbon nanotube - polymer composite using random walks 
###### *Developed by Timothy Burt in the Mullen theoretical condensed matter physics research group at the University of Oklahoma (OU)*
                  
Compatible with Python 2.x/3.x

For serial code, options are given by setting up a config.ini in your choice of save directory
or else default.ini will be read. For parallel code, options are given on the command line.

Serial code is ran by running "python run.py", for the MPI code
run "mpirun -np X python mpi_run.py" with X the number of cores available. **80** cores seems good, or requesting 
10 nodes **exclusively** since a standard node has 8 cores.

Plotting analysis after simulation can be found in plots.py, with stand-alone definitions for individual plots. 
**ipython** can be used to run them.

On OU's Schooner, add this line to your .profile:
`module load OpenMPI mpi4py numpy matplotlib`. `#SBATCH -p 2650` gives a 10 core node, `#SBATCH -p 2670`
gives a 12 core node

Remember you must load conduction before running it, i.e. do not run it in the package directory.


##### References

* Computational modeling of the thermal conductivity of single-walled carbon nanotube-polymer composites
  * H. M. Duong, D. V Papavassiliou, K. J. Mullen, and S. Maruyama, Nanotechnology 19, 065702 (2008).
* Random walks in nanotube composites: Improved algorithms and the role of thermal boundary resistance
  * H. M. Duong, D. V. Papavassiliou, L. L. Lee, and K. J. Mullen, Appl. Phys. Lett. 87, (2005).