#!/bin/bash

cp Nlist24.out clean.out
sed -i 's/PLUMED://g' clean.out
sed -i '/\s\s\s\s\s0/!d' clean.out
sed -i 's/\s\s\s\s\s//g' clean.out
