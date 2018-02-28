#!/bin/bash

if  [ -f "out.out" ] ; then
	rm out.out
fi

gfortran -o potential.o potential.f90
echo 12 | ./potential.o >> out.out

gnuplot --persist -e 'plot "out.out" using 1:3 with line'
gnuplot --persist -e 'plot "out.out" using 1:5 with line'
