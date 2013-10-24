#!/bin/bash
#
# do the grompp 
#
GROMPP="grompp_d"
MDRUN="mdrun_d"
$GROMPP -f md.mdp -c 2ala_ts.gro -p gromacs.top  
$MDRUN -plumed plumed 
