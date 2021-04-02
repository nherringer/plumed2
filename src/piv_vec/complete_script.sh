#!/bin/bash

./run_plumed.sh >& Nlist24.out
./clean.sh
python check_clean.py


