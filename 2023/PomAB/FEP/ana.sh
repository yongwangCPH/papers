#!/bin/sh
#

gmx=/data/public/software/GMX20215/bin/gmx

for i in `seq 0 1 12`
do

cd LAMBDA${i}

echo 0 | $gmx trjconv -f md${i}.gro -o md${i}.pdb -s md${i}.tpr -pbc mol

cd ..
done



