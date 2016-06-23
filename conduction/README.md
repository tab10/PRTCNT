# CONDUCTION
#### Package to study nanotubes in a medium

Program implements a random walk in a medium to calculate its conductivity.
Options are set and specified in default.ini, and the program should be called with "python run.py".

Options are implemented by setting up a config.ini in your choice of save directory
or else default.ini will be read.

Serial code is ran by running "python run.py", for the MPI implementation
run "mpiexec -n X python mpi_run.py" with X the number of cores available.

On OU's Schooner, add this line to your .profile:
`module load OpenMPI mpi4py numpy matplotlib`

##### References

* Computational modeling of the thermal conductivity of single-walled carbon nanotube-polymer composites
  * H. M. Duong, D. V Papavassiliou, K. J. Mullen, and S. Maruyama, Nanotechnology 19, 065702 (2008).
* Random walks in nanotube composites: Improved algorithms and the role of thermal boundary resistance
  * H. M. Duong, D. V. Papavassiliou, L. L. Lee, and K. J. Mullen, Appl. Phys. Lett. 87, (2005).
* Transport of a passive scalar in a turbulent channel flow
  * D. V. D. V. Papavassiliou and T. J. T. J. Hanratty, Int. J. Heat Mass Transf. 40, 1303 (1997).