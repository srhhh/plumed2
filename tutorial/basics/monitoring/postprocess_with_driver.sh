#!/bin/bash
TRJCONV="trjconv_mpi"
DRIVER="plumed driver" 
mv COLVAR COLVAR.bak
#
# generate the pdb
#
echo 0 | $TRJCONV -f 2ala.gro -o 2ala.pdb
echo 0 | $TRJCONV -f traj.trr -o traj.gro
#
# run the driver
#
$DRIVER --pdb 2ala.pdb --igro traj.gro --plumed plumed2.dat 
