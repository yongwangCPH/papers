#!/bin/sh
#

gmx=/data/public/software/GMX20215/bin/gmx
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/public/software/gcc/9.5.0/lib:/data/public/software/gcc/9.5.0/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-11.3/lib:/usr/local/cuda-11.3/lib64

mdp=md.mdp
gro=../gromacs_202208/step7.gro
top=../gromacs_202208/topol.top

#$gmx grompp -f $mdp -c $gro -p $top -o md.tpr  -maxwarn 2 -r $gro -pp all.top
ln -s $gmx PomAB-new
./PomAB-new mdrun -deffnm md -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -gpu_id 2 -cpi -v #-nsteps -1 -update gpu 
