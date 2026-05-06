set terminal pngcairo size 1200,1100 enhanced font 'Arial,11'
set output 'oscillations.png'
set multiplot layout 3,3 title 'Neutrino oscillation probabilities (Earth)'
set ylabel 'E [GeV]'
set xlabel 'cos(zenith)'
set cbrange [0:1]
set palette rgbformulae 33,13,10
set yrange [2:17]
set xrange [-1:0]
set title 'P_{ee}'
plot 'oscillations.dat' using 2:1:3 with image notitle
set title 'P_{e#mu}'
plot 'oscillations.dat' using 2:1:4 with image notitle
set title 'P_{e#tau}'
plot 'oscillations.dat' using 2:1:5 with image notitle
set title 'P_{#mue}'
plot 'oscillations.dat' using 2:1:6 with image notitle
set title 'P_{#mu#mu}'
plot 'oscillations.dat' using 2:1:7 with image notitle
set title 'P_{#mu#tau}'
plot 'oscillations.dat' using 2:1:8 with image notitle
set title 'P_{#taue}'
plot 'oscillations.dat' using 2:1:9 with image notitle
set title 'P_{#tau#mu}'
plot 'oscillations.dat' using 2:1:10 with image notitle
set title 'P_{#tau#tau}'
plot 'oscillations.dat' using 2:1:11 with image notitle
unset multiplot
