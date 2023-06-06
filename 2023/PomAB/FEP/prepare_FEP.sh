#!/bin/sh
#

gmx=/data/public/software/GMX20215/bin/gmx

name=N2K
gro=${name}.gro
# use the conformation that the 1st and 3rd ions have been released
cp ../../gromacs_202208/step7.gro $gro
ln -s ../../gromacs_202208/toppar/ .

top=topol_${name}.top

GPUID=0
for i in `seq 0 1 12`
do

rm -rf LAMBDA${i}
mkdir LAMBDA${i}
cd LAMBDA${i}

cp ../md_FEP.mdp mdrun.mdp
perl -pi -e"s/LAMBDA/$i/g" *.mdp

$gmx grompp -f mdrun.mdp -o md${i}.tpr -c ../$gro -p ../$top -r ../$gro -maxwarn 3 -pp all.top
rm -f \#*

cp ../submit_template.sh submit.sh
perl -pi -e"s/NAME/${name}-${i}/g" submit.sh
perl -pi -e"s/MD/md${i}/g" submit.sh
perl -pi -e"s/GPUID/$GPUID/g" submit.sh

nohup sh submit.sh &

cd ..
done
