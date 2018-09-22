import numpy as np

f = lambda x: ( x - 0.5 )**2

phi = 0.5 * ( 1 + np.sqrt(5) )

a = -1.0
b = 1.0

delta = 0.5 * ( b - a )

eps = 1e-6
while ( delta > eps ):

	p = ( b - a ) / phi
	x1 = b - p
	x2 = a + p

	if ( f(x1) > f(x2) ):
		a = x1
	else:
		b = x2
	delta = 0.5 * ( b - a )

	print('{0:g}\t{1:g}\t{2:g}\t{3:g}'.format( x1, f(x1), x2, f(x2)) )

print( 0.5 * ( a + b ) )
