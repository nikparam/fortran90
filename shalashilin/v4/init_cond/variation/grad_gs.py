import numpy as np
import fort_subroutines as fs
from pprint import pprint
from functools import partial, reduce

def compose(*func):
	init, *rest = reversed(func)
	return lambda *args, **kws: reduce( lambda a, x: x(a), rest, init(*args, **kws) )

lmap = compose( list, map )

switch = 13
params = fs.read_params( switch, '../../potential/params.txt' )


A = 0.88622693
NumG = 4

D = [ -0.0 - 0.517722 * 1j, \
       0.0 + 0.517722 * 1j, \
       0.517722 + 0.0 * 1j, \
      -0.517722 - 0.0 * 1j  ]


omega = [1]*NumG
phase = lmap( lambda x: 0.25 * np.log( x / np.pi ), omega )

a = 0.1 * A
b = 2 * A 
phi = 0.5 * ( 1 + np.sqrt(5) )

delta = abs( b - a )

print(A)

eps = 1.0e-8
while ( delta > eps ):

	p = ( b - a ) / phi
	x1 = b - p
	x2 = a + p

	q1 = [ 0.0, 0.0, -x1, x1 ]
	p1 = [ -x1, x1, 0.0, 0.0 ]

	xi1, eta1 = fs.change_var( q1, p1, omega, phase )
	overlap1 = fs.overlap( xi1, eta1, omega )
	hamiltonian1, dummy = fs.hamiltonian( xi1, eta1, omega, params, overlap1 )
	
	E1 =  np.dot( np.conj(D), np.dot( hamiltonian1, D) )

	q2 = [ 0.0, 0.0, -x2, x2 ]
	p2 = [ -x2, x2, 0.0, 0.0 ]

	xi2, eta2 = fs.change_var( q2, p2, omega, phase )
	overlap2 = fs.overlap( xi2, eta2, omega )
	hamiltonian2, dummy = fs.hamiltonian( xi2, eta2, omega, params, overlap2 )
	
	E2 =  np.dot( np.conj(D), np.dot( hamiltonian2, D) )

	print(x1, E1.real, x2, E2.real)

	if ( E1 > E2 ):
		a = x1
	else:
		b = x2

	delta = abs( b - a )

#	print( '{0:g}\t{1:g}\t{2:g}\t{3:g}'.format(E.real, delta, dE / step, step) )

print()
print( 0.5 * (a+b) )
