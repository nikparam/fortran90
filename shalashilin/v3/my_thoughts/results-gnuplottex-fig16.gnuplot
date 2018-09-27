set terminal latex
set output 'results-gnuplottex-fig16.tex'
set terminal epslatex color size 7, 2.4
set multiplot layout 1, 2

set xlabel 't, a.u.'
set key bottom right
set key at 30., .993
plot "./noneq\_freq/norm\_energy\_1.out" u 1:2 w l ti '$\myint{\Psi}{\Psi}$' lt 1 lw 1 lc rgb 'blue'

set xlabel 't, a.u.'
set ylabel 'E, a.u.'
set yrange [1.268:1.2708]
set key bottom right
set keyt at 30., 1.2685
plot "./noneq\_freq/norm\_energy\_1.out" u 1:3 w l ti '$\mymean{\Psi}{\hat{H}}{\Psi}$' lt 1 lw 1 lc rgb 'red'
