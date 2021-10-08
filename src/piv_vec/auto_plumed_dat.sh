#!/bin/bash

# set PIVATOMS and modify
PIVATOMS=13
sed -i '' "s/PIVATOMS=.*/PIVATOMS=${PIVATOMS}/g" plumed.dat

# creates a string of alternating carbons and a terminating oxygen
# for ATOMTYPES input
string="C1"
for i in $( seq 3 2 23 )
do
	string="${string},C${i}"
done
string="${string},OW"
sed -i '' "s/ATOMTYPES=.*/ATOMTYPES=$string/g" plumed.dat


# set the number of blocks in the PIV
blocks=$(( PIVATOMS*(PIVATOMS-1)/2 ))

# create an array of 1.0 of size blocks
# here it is safely used for both sfactor and sort but independent assignment
# may be required in some cases
declare -a SFACTOR=( $(for i in $( seq 1 $blocks ); do echo 1.0; done) )
sfac=""
for i in "${SFACTOR[@]}"
do
    sfac="$i,$sfac"
done
sfac="${sfac%,}"

sed -i '' "s/SFACTOR=.*/SFACTOR=${sfac}/g" plumed.dat
sed -i '' "s/SORT=.*/SORT=${sfac}/g" plumed.dat

# set NL_CUTOFF and choose a larger radius for C-C than for C-O
for i in {1..25}
do
	ip1=$(( i+1 ))
	for j in $( seq ${ip1} 25)
	do
		if [ $j -eq 25 ]
		then
			val="1.0"
		else
			val="3.5"
		fi
		if [ $i -eq 1 -a $j -eq 1 ]
		then
			NL_CUT="${val}"
		else
			NL_CUT="${NL_CUT},${val}"
		fi
	done
done

sed -i '' "s/NL_CUTOFF=.*/NL_CUTOFF=${NL_CUT}/g" plumed.dat

# set NL_STRIDE and NL_skin in the same way
declare -a NL_STRIDE=( $(for i in $( seq 1 $blocks ); do echo 10.0; done) )
nl_str=""
for i in "${NL_STRIDE[@]}"
do
    nl_str="$i,$nl_str"
done
nl_str="${nl_str%,}"

sed -i '' "s/NL_STRIDE=.*/NL_STRIDE=${nl_str}/g" plumed.dat

declare -a NL_SKIN=( $(for i in $( seq 1 $blocks ); do echo 0.1; done) )
nl_sk=""
for i in "${NL_SKIN[@]}"
do
    nl_sk="$i,$nl_sk"
done
nl_sk="${nl_sk%,}"

sed -i '' "s/NL_SKIN=.*/NL_SKIN=${nl_sk}/g" plumed.dat
