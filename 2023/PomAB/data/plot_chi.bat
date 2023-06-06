set encoding iso_8859_1
set xtics nomirror
set term pdf enhanced color dash size 4.0in, 2.5in

set key out
set xlabel "MD Time (ns)"
set ylabel "Sidechain Rotation (degree)"
set ytics 60

set title "PomA_1"
set output "chi_PomA_1.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainA.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.3 ti "T158-{/Symbol c}_1",\
 '' u 1:($3*180/3.1415926) w lp pt 7 ps 0.3 ti "T185-{/Symbol c}_1",\
 '' u 1:($4*180/3.1415926) w lp pt 7 ps 0.3 ti "T186-{/Symbol c}_1"

set title "PomA_2"
set output "chi_PomA_2.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainB.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.3 ti "T158-{/Symbol c}_1",\
 '' u 1:($3*180/3.1415926) w lp pt 7 ps 0.3 ti "T185-{/Symbol c}_1",\
 '' u 1:($4*180/3.1415926) w lp pt 7 ps 0.3 ti "T186-{/Symbol c}_1"

set title "PomA_3"
set output "chi_PomA_3.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainC.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.3 ti "T158-{/Symbol c}_1",\
 '' u 1:($3*180/3.1415926) w lp pt 7 ps 0.3 ti "T185-{/Symbol c}_1",\
 '' u 1:($4*180/3.1415926) w lp pt 7 ps 0.3 ti "T186-{/Symbol c}_1"

set title "PomA_4"
set output "chi_PomA_4.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainD.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.3 ti "T158-{/Symbol c}_1",\
 '' u 1:($3*180/3.1415926) w lp pt 7 ps 0.3 ti "T185-{/Symbol c}_1",\
 '' u 1:($4*180/3.1415926) w lp pt 7 ps 0.3 ti "T186-{/Symbol c}_1"

set title "PomA_5"
set output "chi_PomA_5.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainE.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.3 ti "T158-{/Symbol c}_1",\
 '' u 1:($3*180/3.1415926) w lp pt 7 ps 0.3 ti "T185-{/Symbol c}_1",\
 '' u 1:($4*180/3.1415926) w lp pt 7 ps 0.3 ti "T186-{/Symbol c}_1"

set title "PomB_1"
set output "chi_PomB_1.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainG.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.3 ti "T21-{/Symbol c}_1",\
 '' u 1:($4*180/3.1415926) w lp pt 7 ps 0.3 ti "D24-{/Symbol c}_2"

# '' u 1:($3*180/3.1415926) w lp pt 7 ps 0.3 ti "D24-{/Symbol c}_1",\

set title "PomB_2"
set output "chi_PomB_2.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainF.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.3 ti "T21-{/Symbol c}_1",\
 '' u 1:($4*180/3.1415926) w lp pt 7 ps 0.3 ti "D24-{/Symbol c}_2"

# '' u 1:($3*180/3.1415926) w lp pt 7 ps 0.3 ti "D24-{/Symbol c}_1",\

set title "Sidechain Rotation of PomA T158"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_T158_2vs5.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainB.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "T158_2",\
 'COLVAR_chainE.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink" ti "T158_5"

set title "Sidechain Rotation of PomA T185"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_T185_2vs5.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainB.txt' u 1:($3*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "T185_2",\
 'COLVAR_chainE.txt' u 1:($3*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink" ti "T185_5"

set title "Sidechain Rotation of PomA T186"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_T186_2vs5.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainB.txt' u 1:($4*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "T186_2",\
 'COLVAR_chainE.txt' u 1:($4*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink" ti "T186_5"

set title "Sidechain Rotation of PomA T33"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_T33.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainA.txt' u 1:($5*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "purple" ti "T33_1",\
 'COLVAR_chainB.txt' u 1:($5*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "T33_2",\
 'COLVAR_chainC.txt' u 1:($5*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "yellow" ti "T33_3",\
 'COLVAR_chainD.txt' u 1:($5*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "green" ti "T33_4",\
 'COLVAR_chainE.txt' u 1:($5*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink" ti "T33_5"

set title "Sidechain Rotation of PomA T158"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_T158.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainA.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "purple" ti "T158_1",\
 'COLVAR_chainB.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "T158_2",\
 'COLVAR_chainC.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "yellow" ti "T158_3",\
 'COLVAR_chainD.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "green"  ti "T158_4",\
 'COLVAR_chainE.txt' u 1:($2*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink"   ti "T158_5"

set title "Sidechain Rotation of PomA T185"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_T185.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainA.txt' u 1:($3*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "purple" ti "T185_1",\
 'COLVAR_chainB.txt' u 1:($3*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "T185_2",\
 'COLVAR_chainC.txt' u 1:($3*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "yellow" ti "T185_3",\
 'COLVAR_chainD.txt' u 1:($3*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "green"  ti "T185_4",\
 'COLVAR_chainE.txt' u 1:($3*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink"   ti "T185_5"

set title "Sidechain Rotation of PomA T33"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_T33_2vs5.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainB.txt' u 1:($5*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "T33_2",\
 'COLVAR_chainE.txt' u 1:($5*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink" ti "T33_5"

set title "Sidechain Rotation of PomA M155"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_M155.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainA.txt' u 1:($6*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "purple" ti "M155_1",\
 'COLVAR_chainB.txt' u 1:($6*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "M155_2",\
 'COLVAR_chainC.txt' u 1:($6*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "yellow" ti "M155_3",\
 'COLVAR_chainD.txt' u 1:($6*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "green" ti "M155_4",\
 'COLVAR_chainE.txt' u 1:($6*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink" ti "M155_5"

set title "Sidechain Rotation of PomA M155"
set ylabel "{/Symbol c}_1 (degree)"
set output "chi_PomA_M155_2vs5.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainB.txt' u 1:($6*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "orange" ti "M155_2",\
 'COLVAR_chainE.txt' u 1:($6*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "pink" ti "M155_5"

set title "Dynamics of PomB D24"
set ylabel "{/Symbol c}_2 (degree)"
set output "chi_PomB_D24_1vs2.pdf"
p [0:1000][-180:180] \
 'COLVAR_chainG.txt' u 1:($4*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "grey" ti "D24_1",\
 'COLVAR_chainF.txt' u 1:($4*180/3.1415926) w lp pt 7 ps 0.2 lc rgb "black" ti "D24_2"


