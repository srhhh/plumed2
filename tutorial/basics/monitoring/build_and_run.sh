#!/bin/bash
#
# do the grompp 
#
grompp_mpi -f md.mdp -c 2ala.gro -p gromacs.top  
mdrun_mpi -plumed plumed 
