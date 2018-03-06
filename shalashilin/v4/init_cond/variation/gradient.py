import numpy as np
import fort_subroutines as fs
from pprint import pprint
from functools import partial, reduce
import matplotlib.pyplot as plt

def compose(*func):
	init, *rest = reversed(func)
	return lambda *args, **kws: reduce( lambda a, x: x(a), rest, init(*args, **kws) )

def normalization( q, p, D, omega, phase ):
	xi, eta = fs.change_var( q, p, omega, phase )
	overlap = fs.overlap( xi, eta, omega )
	N = np.dot( np.conj(D), np.dot( overlap, D ) )
	D /= np.sqrt( N )

	return D

def num_zeroes( func, x_array ):

	a = x_array[0]
	NPoints = len( x_array )
	step = 2 * a / NPoints

	df = partial( lambda f, x1, dx: ( f( x1 + dx ) - f( x1 ) ) / dx, func )

	df0 = 0.0
	df1 = 0.0
	count = 0
	for j in x_array:
		df0 = df(j,step)
		if df1 < 0 and df0 > 0: 
			count += 1
		df1 = df0

	return count


lmap = compose( list, map )

switch = 13
params = fs.read_params( switch, '../../potential/params.txt' )


A = 0.9
NumG = 8

#D = [ 4.86761400 - 12.14180620 * 1j, \
#      14.45370904 - 31.39053856 * 1j, \
#      -0.03538665 + 26.55952747 * 1j, \
#      -19.28601170 + 16.97285778 * 1j ]

#D = [  23.29077405 - 104.08655773 * 1j, \
#      -23.00668126 + 104.13895337 * 1j, \
#       89.79251688 + 30.71022183 * 1j, \
#      -90.07539352 - 30.76206115 * 1j ]

#D = [ 1, 1, 1, 1, 1j, 1j, 1j, 1j ]
#D = [ -1j, 1j, 1, -1 ]
#D = [ 1, 1, -1, -1 ]
#D = [ -1j, 1j, -1, 1 ]

q = [ 0.0, 0.0, -A, A, -0.5 * A, -0.5 * A, 0.5 * A, 0.5 * A ]
p = [ -A, A, 0.0, 0.0, -0.5 * A, 0.5 * A, -0.5 * A, 0.5 * A ]

omega = [1]*NumG
phase = lmap( lambda x: 0.25 * np.log( x / np.pi ), omega )

D = normalization( q, p, D, omega, phase )

psi = lambda c, q0, p0, omega, phase, x: c * np.exp( - 0.5 * omega * ( x -q0 )**2 + \
					 1j * p0 * ( x - q0 ) + phase )
mapped_psi = lmap( lambda c, q0, p0, o, p: partial( psi, c, q0, p0, o, p ), \
		   					 D, q, p, omega, phase )
full_psi = lambda x: reduce( lambda a, f: a + f(x), mapped_psi[1:], mapped_psi[0](x) )
abs_full_psi = lambda x: abs( full_psi(x) )**2

a = 5
NPoints = 1000
x = np.linspace(-a,a,NPoints)

NumZ = num_zeroes( abs_full_psi, x )

fig = plt.figure( figsize = ( 9, 9 ), facecolor = 'w' )
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(224)

ax1.plot( x, lmap( abs_full_psi, x ) )

eps = 1.0e-8

E_exact = 2.5
E0 = 2 * eps
delta = 2 * eps
step = -0.01

der_E = 2 * eps

while ( abs( delta ) > eps ):
	xi, eta = fs.change_var( q, p, omega, phase )
	overlap = fs.overlap( xi, eta, omega )
	D = normalization( q, p, D, omega, phase )

	hamiltonian, dummy = fs.hamiltonian( xi, eta, omega, params, overlap )
	E =  np.dot( np.conj(D), np.dot( hamiltonian, D) )

	delta = E.real - E0.real
	delta_a = E.real - E_exact
	der_E = delta / step

	print( 'zeroes = {0:g}\tA = {1:g}\tE = {2:g}\tstep = {3:g}\tE-E_ex = {4:g}\tdE = {5:g}\tE\' = {6:g}'\
		.format( NumZ, A, E.real, step, delta_a, delta, der_E ) )

	tempNumZ = num_zeroes( abs_full_psi, x )
	if ( tempNumZ != NumZ ):
		step *= -0.1 
	else:
		step *= 1 

	E0 = E
	A += step

	q = [ 0.0, 0.0, -A, A, -0.5 * A, -0.5 * A, 0.5 * A, 0.5 * A ]
	p = [ -A, A, 0.0, 0.0, -0.5 * A, 0.5 * A, -0.5 * A, 0.5 * A ]
	mapped_psi = lmap( lambda c, q0, p0, o, p: partial( psi, c, q0, p0, o, p ), \
			   					 D, q, p, omega, phase )
	full_psi = lambda x: reduce( lambda a, f: a + f(x), mapped_psi[1:], mapped_psi[0](x) )
	abs_full_psi = lambda x: abs( full_psi(x) )**2


ax1.plot( x, lmap( abs_full_psi, x ) )
ax1.text( 0, 0.7, r'$|\psi|^2$', fontsize = 20 )

for i in range( NumG ):
	ax2.plot( x, lmap( lambda x: x.real, lmap( mapped_psi[i], x ) ) )
	ax3.plot( x, lmap( lambda x: x.imag, lmap( mapped_psi[i], x ) ) )

ax2.text( -0.8, -0.05, r'$\Re\psi$', fontsize = 20 )
ax3.text( -0.8, -0.035, r'$\Im\psi$', fontsize = 20 )

print('\n{0:g}'.format(A))
print(D)
print(q,p)
plt.show()