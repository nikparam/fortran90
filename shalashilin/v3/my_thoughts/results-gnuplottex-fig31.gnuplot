set terminal latex
set output 'results-gnuplottex-fig31.tex'
set terminal epslatex color size 7, 2.4
set multiplot layout 1, 2

set xlabel 't, a.u.'
set ylabel 'E, a.u.'
set key bottom right
set keyt at 16., .995
plot "./exp/norm\_energy.out" u 1:3 w l ti '$\mymean{\Psi}{\hat{H}}{\Psi}$' lt 1 lw 1 lc rgb 'red'

set xlabel 't, a.u.'
set ylabel 'E, a.u.'
set ytics 0.01
set key bottom right
set keyt at 16., 2.235
plot "./exp/norm\_energy.out" u 1:8 w l ti '$\mathit{Tr}(\mathbbm{H})$' lt 1 lw 1 lc rgb 'blue'