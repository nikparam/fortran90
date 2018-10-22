set terminal latex
set output 'results-gnuplottex-fig33.tex'
set terminal epslatex color size 7, 2.4
set multiplot layout 1, 2

set xlabel '<q>, a.u.'
set ylabel '$\frac{\partial \text{V}}{\partial \text{x}}$, a.u.'
set key bottom right
set ytics 0.2
set yrange [-1.2:0.1]
set keyt at 20., -1.
plot "./trace/main.log" u 2:3 w l ti '$\underset{k=1}{\overset{2}{\sum}}\mymean{g_k}{V(x-q_k)}{g_k}$' lt 1 lw 1 lc rgb 'red'

set xlabel 't, a.u.'
set ylabel '$\frac{\partial \text{V}}{\partial \text{x}}$, a.u.'
set key spacing 1.75
set key bottom right
set ytics 0.2
set keyt at 16., -1.2
set yrange [-1.4:0.1]
plot "./trace/main.log" u 1:3 w l ti '$\underset{k=1}{\overset{2}{\sum}}\mymean{g_k}{V(x-q_k)}{g_k}$' lt 1 lw 1 lc rgb 'red'
