set terminal pngcairo size 1200,1100 enhanced font 'Arial,11'
set output 'oscillations.png'
set multiplot layout 3,3 title 'Neutrino oscillation probabilities (Earth)'
set xlabel 'E [GeV]'
set ylabel 'cos(zenith)'
set cbrange [0:1]
set palette rgbformulae 33,13,10
set xrange [2:40]
set yrange [-1:0]
set title 'P_{ee}'
plot 'oscillations.dat' using 1:2:3 with image notitle
set title 'P_{e#mu}'
plot 'oscillations.dat' using 1:2:4 with image notitle
set title 'P_{e#tau}'
plot 'oscillations.dat' using 1:2:5 with image notitle
set title 'P_{#mue}'
plot 'oscillations.dat' using 1:2:6 with image notitle
set title 'P_{#mu#mu}'
plot 'oscillations.dat' using 1:2:7 with image notitle
set title 'P_{#mu#tau}'
plot 'oscillations.dat' using 1:2:8 with image notitle
set title 'P_{#taue}'
plot 'oscillations.dat' using 1:2:9 with image notitle
set title 'P_{#tau#mu}'
plot 'oscillations.dat' using 1:2:10 with image notitle
set title 'P_{#tau#tau}'
plot 'oscillations.dat' using 1:2:11 with image notitle
unset multiplot
