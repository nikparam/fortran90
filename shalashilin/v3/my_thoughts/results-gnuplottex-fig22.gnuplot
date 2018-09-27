set terminal latex
set output 'results-gnuplottex-fig22.tex'
set terminal epslatex color

set size 1, 0.7
set xlabel 't, a.u.'
set ylabel 'E, a.u.'
set key bottom right
set yrange [-1.75:1.75]
set keyt at 100., -1.4
plot "./noneq\_freq/norm\_energy\_2.out" u 1:4 w l ti '$\mymean{\Psi}{\hat{T} - \hat{V}}{\Psi}$' lt 1 lw 1 lc rgb 'green'
