set terminal latex
set output 'results-gnuplottex-fig21.tex'
set terminal epslatex color size 7, 2.4
set multiplot layout 1, 2

set xlabel 't, a.u.'
set key bottom right
set key at 100., .993
plot "./noneq\_freq/norm\_energy\_2.out" u 1:2 w l ti '$\myint{\Psi}{\Psi}$' lt 1 lw 1 lc rgb 'blue'

set xlabel 't, a.u.'
set ylabel 'E, a.u.'
set key bottom right
set keyt at 100., 1.26
plot "./noneq\_freq/norm\_energy\_2.out" u 1:3 w l ti '$\mymean{\Psi}{\hat{H}}{\Psi}$' lt 1 lw 1 lc rgb 'red'
