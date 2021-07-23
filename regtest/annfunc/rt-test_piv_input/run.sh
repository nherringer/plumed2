#!/bin/bash

#plumed driver --plumed plumed.dat --mf_pdb md.pdb >& out.out

export OMP_NUM_THREADS=5
#gmx mdrun -ntomp $OMP_NUM_THREADS -s ../../../../codes/simple_run_template/3_md/md.tpr -deffnm plumed-test-md -pin on -plumed plumed.dat -v

plumedFile=plumed-ann-full.dat
gmx mdrun -ntomp $OMP_NUM_THREADS -s ../../../../codes/simple_run_template/3_md/md.tpr -deffnm plumed-test-md -pin on -plumed $plumedFile -v

