#!/bin/sh
#

gmx=/data/public/software/GMX20215/bin/gmx

#echo 1 18 | $gmx trjconv -f md.cpt -o align.pdb -s md.tpr -fit rot+trans -n index.ndx

/data/public/software/GROMAPs/bin/gmx maptide -f align.xtc -s align.pdb -select 'resname SOD' -spacing 0.1 -mo SOD_dens.ccp4 

/data/public/software/GROMAPs/bin/gmx maptide -f align.xtc -s align.pdb -select 'resname TIP3' -spacing 0.1 -mo W_dens.ccp4 

rm -f \#*
