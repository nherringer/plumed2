#!/bin/bash

#plumed driver --plumed plumed.dat --mf_pdb md.pdb >& out.out

export OMP_NUM_THREADS=5
#gmx mdrun -ntomp $OMP_NUM_THREADS -s ../../../../codes/simple_run_template/3_md/md.tpr -deffnm plumed-test-md -pin on -plumed plumed.dat -v

#cd ../../../../../codes/simple_run_template/
#source run_sim.sh
#cd /data/research/enhancements/plumed-development/plumed2/regtest/annfunc/rt-test_piv_input/sys_4_at_9_piv_updated

plumedFile=plumed.dat
mpirun -np 3 gmx mdrun -ntomp $OMP_NUM_THREADS -s ../../../../../codes/simple_run_template/3_md/md.tpr -pin on -deffnm plumed-test-md -plumed $plumedFile 
