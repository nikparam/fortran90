import numpy as np
from itertools import chain
from functools import reduce
from pprint import pprint

def compose(*func):
	init, *rest = reversed(func)
	return lambda *args, **kws: reduce( lambda a, x: x(a), rest, init(*args, **kws) )

def string_out(seq):
	return reduce( lambda a, x: a + ' & ' + x , seq[1:], seq[0] ) + '\\\ \n'

floatify = lambda x: float(x)
tostr = lambda array: list( map( lambda x: str( round(x,4) ), array ) )
f = lambda x, y: x**2 + y**2
g = lambda a, x: a + x

lmap = compose( list, map )

with open( '../eigen_vectors.out', 'r' ) as fin:
	NumStates = int( fin.readline().split()[0] )
	energy = list( map( floatify, fin.readline().split() ) )
	eigen_vectors = []
	for _ in range( NumStates ):
		vector = np.transpose( np.reshape( lmap( floatify, fin.readline().split() ), (NumStates, 2) ) )
		eigen_vectors.append( lmap( f, vector[0], vector[1] ) )		

eigen_vectors = np.transpose( eigen_vectors )

with open( '../qp_init.out', 'r' ) as fin:
	q = []
	p = []
	for _ in range(3):
		fin.readline()
	for _ in range(NumStates):
			line = fin.readline()
			q.append( float( line.split()[0] ) )
			p.append( float( line.split()[1] ) )

with open( 'weights.out', 'w' ) as fout:
	norm_eigen_vectors = []
	for _ in range( NumStates ):
		norm_eigen_vectors.append(  lmap( lambda x: x / reduce( g, eigen_vectors[_] ), eigen_vectors[_] ) )

	norm_eigen_vectors = np.transpose( norm_eigen_vectors )
	temp = ['q','p']
	temp.extend( tostr( energy ) )
	fout.write( string_out( temp ) ) 
	for _ in range(NumStates):
		temp = [ q[_], p[_] ]
		temp.extend( norm_eigen_vectors[_][:] )
		fout.write( string_out( tostr( temp ) ) )


