set style line 1 lt 1 lw 1 lc rgb 'red'
set style line 2 lt 1 lw 1 lc rgb 'blue'

set style arrow 1 nohead ls 1
set style arrow 2 nohead ls 2

#########################################################################################
# plot grid
#########################################################################################

set term x11 0
set xrange [-5 : 5]
set autoscale y
plot '../qp_init.out' every ::3 using 1:2

#########################################################################################
# plot deviation of energy as function of calculated energy
#########################################################################################

set term x11 1

set autoscale x

plot '../energy_states.out' using 1:4

#########################################################################################
# plot energy levels --- red for Harmonic oscillator, blue for calc. energies
#########################################################################################

set term x11 2

set xtics format " "

plot '../energy_states.out' using (0.0):1:(1.0):(0.0) with vectors arrowstyle 1, \
     '../energy_states.out' using (0.0):2:(1.0):(0.0) with vectors arrowstyle 2


