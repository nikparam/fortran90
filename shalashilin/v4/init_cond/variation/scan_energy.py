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


A = np.linspace( 0.2, 2.0, 25 )
NumG = 4

D = [ 1.0, 1.0, 1.0, 1.0 ]

omega = [1]*NumG
phase = lmap( lambda x: 0.25 * np.log( x / np.pi ), omega )

for a in A:
	q = [ 0.0, 0.0, -a, a ]
	p = [ -a, a, 0.0, 0.0 ]

	xi, eta = fs.change_var( q, p, omega, phase )
	overlap = fs.overlap( xi, eta, omega )

	N =  np.dot( np.conj(D), np.dot( overlap, D) )
	D /= np.sqrt(N)

	hamiltonian, dummy = fs.hamiltonian( xi, eta, omega, params, overlap )
	
	E =  np.dot( np.conj(D), np.dot( hamiltonian, D) )

	print( a, E.real )

