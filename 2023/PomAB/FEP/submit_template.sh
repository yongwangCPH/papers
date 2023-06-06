#!/bin/sh
#

gmx=/data/public/software/GMX20215/bin/gmx
export GMX_ALLOW_CPT_MISMATCH=1
export PLUMED_USE_LEPTON=yes
export PLUMED_KERNEL="/data/public/software/PLUMED274/lib/libplumedKernel.so"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/public/software/openmpi4/lib/

run=NAME
ln -s $gmx $run
./$run mdrun -deffnm MD -ntmpi 1 -ntomp 2 -nb gpu -pme gpu -bonded gpu -gpu_id GPUID -cpi -v 
