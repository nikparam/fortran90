gnuplot --persist -e 'plot "energy01.out" using 1:4 with line'
gnuplot --persist -e 'plot "out01.out" using 1:2 with line'
gnuplot --persist -e 'plot "out01.out" using 1:3 with line'
gnuplot --persist -e 'plot "out01.out" using 2:3 with line'

