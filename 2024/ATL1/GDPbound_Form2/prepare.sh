#!/bin/sh
#

gmx=/data/public/software/GMX20225/bin/gmx

gro=hATL_solv_ions.gro 
top=topol.top
ndx=index.ndx

export GMX_MAXCONSTRWARN=-1

pre=1
if [ $pre -eq 1 ]; then
# step6.0
$gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c $gro -r $gro -p $top -pp all.top -maxwarn 2
$gmx mdrun -v -deffnm step6.0_minimization -ntmpi 1 -ntomp 8 -gpu_id 0

# Equilibration
cnt=1
cntmax=6

while [ ${cnt} -le ${cntmax} ]
do
    pcnt=`echo $cnt | awk '{print $1-1}'`

    if [ ${cnt} -eq 1 ]; then
        $gmx grompp -f step6.${cnt}_equilibration.mdp -o step6.${cnt}_equilibration.tpr -c step6.${pcnt}_minimization.gro -r $gro  -p $top -maxwarn -1 -n $ndx
        $gmx mdrun -v -deffnm step6.${cnt}_equilibration -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -gpu_id 0
    else
        $gmx grompp -f step6.${cnt}_equilibration.mdp -o step6.${cnt}_equilibration.tpr -c step6.${pcnt}_equilibration.gro -r $gro  -p $top -maxwarn -1 -n $ndx
        $gmx mdrun -v -deffnm step6.${cnt}_equilibration -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -gpu_id 0
    fi

    cnt=$(( cnt+1 ))

done
fi

# Production
$gmx grompp -f step7_production.mdp -o step7.tpr -c step6.6_equilibration.gro  -p $top -r step6.6_equilibration.gro -n $ndx -maxwarn 2
$gmx mdrun -v -deffnm step7 -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -gpu_id 0
