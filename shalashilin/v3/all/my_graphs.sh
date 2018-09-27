set terminal qt 0
set title '<Psi|Psi>' font ',15'
plot "norm_energy.out" u 1:2 w l

set terminal qt 1
set title '<H>' font ',15'
plot "norm_energy.out" u 1:3 w l

set terminal qt 2
set title '<T> -- <V>' font ',15'
plot "norm_energy.out" u 1:4 w l

set terminal qt 3
set title '0.5 * <p>**2 / m + <V>' font ',15'
plot "norm_energy.out" u 1:5 w l
