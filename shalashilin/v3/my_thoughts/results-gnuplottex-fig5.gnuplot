set terminal latex
set output 'results-gnuplottex-fig5.tex'
set terminal epslatex color

set size 1, 0.6
set xlabel 't, a.u.'
set key bottom right
set yrange [0.48:0.515]
set keyt at 10., .481
plot "./equal\_freq/width\_1.out" u 1:2 w l ti '$\langle x^2\rangle-\langle x\rangle^2$' lt 1 lw 1 lc rgb 'purple'
