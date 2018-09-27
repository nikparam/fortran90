set terminal latex
set output 'results-gnuplottex-fig18.tex'
set terminal epslatex color

set size 1, 0.75
set xlabel 't, a.u.'
set ylabel 'E, a.u.'
set yrange [0.99:1.01]
set key bottom right
set keyt at 30., .995
plot "./noneq\_freq/norm\_energy\_1.out" u 1:5 w l ti '$\frac{1}{2m}\mymean{\Psi}{\hat{p}}{\Psi}^2 + \frac{1}{2}m\omega^2\mymean{\Psi}{\hat{q}}{\Psi}^2$' lt 1 lw 1 lc rgb 'orange'
