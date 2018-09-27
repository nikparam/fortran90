set terminal latex
set output 'results-gnuplottex-fig9.tex'
set terminal epslatex color size 7, 3
set multiplot layout 1, 2

set xlabel 'q, a.u'
set ylabel 'p, a.u.'
set key bottom right
set key at 0.75, -0.125
plot "./equal\_freq/coords\_2.out" u 2:5 w l ti '$p_1 : q_1$' lt 1 lw 2 lc rgb 'blue', \
     "./equal\_freq/coords\_2.out" u 3:6 w l ti '$p_2 : q_2$' lt 1 lw 2 lc rgb 'red', \
     "./equal\_freq/coords\_2.out" u 4:7 w l ti '$p_3 : q_3$' lt 1 lw 2 lc rgb 'green', \

set xlabel 'q, a.u'
set ylabel 'p, a.u.'
set key bottom right
set key at 0.75, -0.125
plot "./equal\_freq/qp\_mean\_2.out" u 2:3 w l ti '$\langle p\rangle : \langle q\rangle$' lt 3 lw 2 lc rgb 'purple'
