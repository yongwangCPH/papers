### 1D FES

# sigma for KDE      
sigma1=0.2206860033617057
sigma2=0.08818475008330701  
sigma3=0.1
sigma4=0.1

# CV name
cv1=tica0
cv2=tica1
cv3=phi
cv4=psi

# COLVAR file suffix
k=0
# FES printing stride := n_sim * intended_stride
stride=4800

#mkdir ${k}_error
# block analysis
for i in 2 4 10 15 30 50 75 100 200 400 600 800 1000
do
    python3 FES_from_Reweighting.py --colvar COLVAR_${k}_sort --outfile fes-$j.dat --sigma $sigma1 --temp 310 --cv $cv1 --block $i
    rm fes*
    python3 FES_from_Reweighting.py --colvar COLVAR_${k}_sort --outfile fes-$j.dat --sigma $sigma2 --temp 310 --cv $cv2 --block $i
    rm fes*
    python3 FES_from_Reweighting.py --colvar COLVAR_${k}_sort --outfile fes-$j.dat --sigma $sigma3 --temp 310 --cv $cv3 --block $i
    rm fes*
    python3 FES_from_Reweighting.py --colvar COLVAR_${k}_sort --outfile fes-$j.dat --sigma $sigma4 --temp 310 --cv $cv4 --block $i
    rm fes*
done

mkdir ${k}_error
mv err* ${k}_error

# Print FES
python3 FES_from_Reweighting.py --colvar COLVAR_${k}_sort --outfile fes-$cv1.dat --sigma $sigma1 --temp 310 --cv $cv1 --stride $stride
python3 FES_from_Reweighting.py --colvar COLVAR_${k}_sort --outfile fes-$cv2.dat --sigma $sigma2 --temp 310 --cv $cv2 --stride $stride
python3 FES_from_Reweighting.py --colvar COLVAR_${k}_sort --outfile fes-$cv3.dat --sigma $sigma3 --temp 310 --cv $cv3 --stride $stride 
python3 FES_from_Reweighting.py --colvar COLVAR_${k}_sort --outfile fes-$cv4.dat --sigma $sigma4 --temp 310 --cv $cv4 --stride $stride

mkdir hist-fes$k
mv fes* hist-fes$k/

mv log$cv1 log${cv1}_${k}
mv log$cv2 log${cv2}_${k}
mv log$cv3 log${cv3}_${k}
mv log$cv4 log${cv4}_${k}
