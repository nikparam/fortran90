gnuplot --persist -e 'plot "out.out" u 1:5 w l, "out.out" u 1:6 w l'
gnuplot --persist -e 'plot "out.out" u 1:7 w l, "out.out" u 1:8 w l, "out.out" u 1:($5+$6) w l'
gnuplot --persist -e 'plot "out.out" u 1:4 w l'


