set terminal latex
set output 'results-gnuplottex-fig30.tex'
set terminal epslatex color

set size 1, 0.6
set xlabel 't, a.u.'
set key bottom right
set yrange [0.2:1.7]
set keyt at 100., .3
plot "./noneq\_freq/width\_3.out" u 1:2 w l ti '$\langle x^2\rangle-\langle x\rangle^2$' lt 1 lw 1 lc rgb 'purple'
