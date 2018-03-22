import numpy as np
import fort_subroutines as fs
from functools import partial, reduce
from itertools import product
from pprint import pprint

def compose(*func):
	init, *rest = reversed(func)
	return lambda *args, **kws: reduce( lambda a, x: x(a), rest, init(*args, **kws) )

lmap = compose( list, map )

def print_matrix( matrix ):
	for row in matrix:
		line = ''
		for element in row:
			line += '{:.4f}'.format( element ) + ' '
		print(line)
	
	print()


def orth_overlap( S, i, j ):
	N = overlap[i][j]
	for k in range(i):
		N -= S[ i ][ k ] * S[ k ][ i ]
	for m in range(j):
		N -= S[ j ][ m ] * S[ m ][ j ]
	for k, m in product( range( i ), range( j ) ):
		N += S[ i ][ k ] * S[ k ][ m ] * S[ m ][ j ]
	return N
			

switch = 13
params = fs.read_params( switch, '../../potential/params.txt' )

A = 0.25
NumG = 4

omega = [1]*NumG
phase = lmap( lambda x: 0.25 * np.log( x / np.pi ), omega )

q = [ 0.0, 0.0, -A, A ]
p = [ -A, A, 0.0, 0.0 ]

xi, eta = fs.change_var( q, p, omega, phase )
overlap = fs.overlap( xi, eta, omega )
hamiltonian, dummy = fs.hamiltonian( xi, eta, omega, params, overlap )

S = [ [0] * NumG ] * NumG
for i in range( NumG ):
	S[i][i] = orth_overlap( overlap, i, i ).real
	for j in range( i ):
		S[i][j] = orth_overlap( overlap, i, j ) / np.sqrt( S[i][i] )
		S[j][i] = np.conj( S[i][j] )

print_matrix( S )


