### 2D FES

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
stride=40000000000000

python3 FES_from_Reweighting.py --colvar COLVAR_$k --outfile fes2d-${cv1}-${cv2}.dat --sigma $sigma1,$sigma2 --temp 310 --cv $cv1,$cv2 --stride $stride
python3 FES_from_Reweighting.py --colvar COLVAR_$k --outfile fes2d-${cv3}-${cv4}.dat --sigma $sigma3,$sigma4 --temp 310 --cv $cv3,$cv4 --stride $stride

mkdir hist-fes$k
mv fes* hist-fes$k/

