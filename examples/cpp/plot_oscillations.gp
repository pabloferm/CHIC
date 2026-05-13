set terminal pngcairo size 1200,1100 enhanced font 'Arial,11'
set output 'oscillations.png'
set multiplot layout 3,3 title 'Neutrino oscillation probabilities (Earth)'
set ylabel 'E [GeV]'
set xlabel 'cos(zenith)'
set palette rgbformulae 33,13,10
set yrange [2.0:17]
set xrange [-1:0]
unset yrange
unset xrange
stats 'oscillations.dat' using 3 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{ee}'
plot 'oscillations.dat' using 2:1:3 with image notitle
unset yrange
unset xrange
stats 'oscillations.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{e#mu}'
plot 'oscillations.dat' using 2:1:4 with image notitle
unset yrange
unset xrange
stats 'oscillations.dat' using 5 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{e#tau}'
plot 'oscillations.dat' using 2:1:5 with image notitle
unset yrange
unset xrange
stats 'oscillations.dat' using 6 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{#mue}'
plot 'oscillations.dat' using 2:1:6 with image notitle
unset yrange
unset xrange
stats 'oscillations.dat' using 7 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{#mu#mu}'
plot 'oscillations.dat' using 2:1:7 with image notitle
unset yrange
unset xrange
stats 'oscillations.dat' using 8 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{#mu#tau}'
plot 'oscillations.dat' using 2:1:8 with image notitle
unset yrange
unset xrange
stats 'oscillations.dat' using 9 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{#taue}'
plot 'oscillations.dat' using 2:1:9 with image notitle
unset yrange
unset xrange
stats 'oscillations.dat' using 10 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{#tau#mu}'
plot 'oscillations.dat' using 2:1:10 with image notitle
unset yrange
unset xrange
stats 'oscillations.dat' using 11 nooutput
set cbrange [STATS_min:STATS_max]
set yrange [2.0:17]
set xrange [-1:0]
set title 'P_{#tau#tau}'
plot 'oscillations.dat' using 2:1:11 with image notitle
