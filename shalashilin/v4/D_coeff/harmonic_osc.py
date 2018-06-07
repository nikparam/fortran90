import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
import scipy.integrate as spint

key = 0

def lagrangian(q,p,m,omega):

	return 0.5 * p**2 / m - 0.5 * m * omega**2 * q**2

def overlap( q, p, m, omega, norm, N ):

	xi = lambda x, y: m * omega * x + 1j * y
	eta = lambda x, y: 0.25 * np.log(norm) - 0.5 * m * omega * x**2 -1j * x * y

	overlap = []

	for i, j in zip( map(xi,q,p), map(eta,q,p) ):
		for k, l in zip( map(xi,q,p), map(eta,q,p) ):
			overlap.append( np.sqrt(1/norm) * np.exp( (np.conjugate(i)+k)**2/(4 * m * omega ) + np.conjugate(j)+l) )	

	return np.reshape(overlap,[N,N])

def kinetic( q, p, S, m, omega, N ):

	xi = lambda x, y: m * omega * x + 1j * y

	kinetic = []

	for i in map(xi,q,p):
		for j in map(xi,q,p):
			kinetic.append( 0.5 * omega - 0.5 * j**2 / m + 0.5 * j * (np.conjugate(i)+j) / m - #
					0.5 * m * omega**2 * ( 0.5 / ( m * omega ) + ( 0.5 * (np.conjugate(i)+j) / ( m * omega) )**2 ) )

	return np.reshape(kinetic,[N,N]) * S

def potential( q, p, S, m, omega, N ):

	xi = lambda x, y: m * omega * x + 1j * y

	potential = []

	for i in map(xi,q,p):
		for j in map(xi,q,p):
			potential.append( 0.5 * m * omega**2 * ( 0.5 / ( m * omega ) + ( 0.5 * (np.conjugate(i)+j) / ( m * omega) )**2 ) )

	return np.reshape(potential,[N,N]) * S

def hamiltonian( q, p, S, m, omega, N ):

	xi = lambda x, y: m * omega * x + 1j * y

	hamiltonian = []

	for i in map(xi,q,p):
		for j in map(xi,q,p):
			hamiltonian.append( 0.5 * ( omega + np.conjugate(i) * j / m ) )

	return np.reshape(hamiltonian,[N,N]) * S

def full_lagrangian( q, p, S, m, omega, N ):

	xi = lambda x, y: m * omega * x + 1j * y

	lagrangian = []

	for i in map(xi,q,p):
		for j in map(xi,q,p):
			lagrangian.append( -0.25 * (np.conjugate(i)**2+j**2) / m )

	return np.reshape(lagrangian,[N,N]) * S

def func( t, y, args ):

	ydot = [0]*3*args[2]

	for i in range( args[2] ):
		ydot[i] = y[ args[2] + i ] / args[0]
		ydot[i+args[2]] = -args[0] * args[1]**2 * y[ i ] 
		ydot[i+2*args[2]] = -1j*( 0.5 * args[1] - lagrangian( y[ i ], y[ args[2] + i ], args[0], args[1] ) )*y[ 2 * args[2] + i ]

	return ydot
	

A = 3.56072
m = 1.0
omega = 1.0
r = 10
n = 2
norm = m * omega / np.pi

q = [A, -A, -A, A, 2*A, 0.0, -2*A, 0.0]
p = [A, A, -A, -A, 0.0, 2*A, 0.0, -2*A]

#q = [A, 0.0, -A, 0.0]
#p = [0.0, A, 0.0, -A]

if len(q) == len(p):
	N = len(q)
else:
	raise NameError('initial conditions have inconsistent length')

C = [1.0]*N

S = overlap(q,p,m,omega,norm,N)
s = np.dot( np.conjugate(C),np.dot(S,C))
C /= np.sqrt(s)

y0 = q[:]
y0.extend(p)
y0.extend(C)

t0 = 0.
res = spint.ode(func).set_integrator('zvode', method='bdf',rtol=1e-11, atol=1e-11)
res.set_initial_value(y0,t0).set_f_params([m,omega,N,r,n])

t1 = 10
dt = 0.05

if key == 0:
	e = np.dot( np.conjugate(C),np.dot(hamiltonian(q,p,S,m,omega,N),C))
	l = np.sqrt(2 * e.real / m / omega**2)
	x = np.linspace(-l-5,l+5,1000)
	ax = plt.gca()
	ax.set_ylim([e.real - 0.5, e.real + 0.5])
	ax.set_xlim( [ -l - 5, l + 5 ] )
	ax.plot(x,[0.5 * m * omega**2 * _**2 for _ in x])
	psi = [0]*len(x)
	current_plot1, = ax.plot( x, [abs(_)**2 for _ in psi] )

while res.successful() and res.t < t1:

	y = res.integrate( res.t + dt )
	q = y[:N]
	p = y[N:2*N]
	C = y[2*N:]

	S = overlap(q,p,m,omega,norm,N)
	s = np.dot( np.conjugate(C),np.dot(S,C))
	e = np.dot( np.conjugate(C),np.dot(hamiltonian(q,p,S,m,omega,N),C))
	l = np.dot( np.conjugate(C),np.dot(full_lagrangian(q,p,S,m,omega,N),C))
	t = np.dot( np.conjugate(C),np.dot(kinetic(q,p,S,m,omega,N),C))
	v = np.dot( np.conjugate(C),np.dot(potential(q,p,S,m,omega,N),C))
	
	if key == 0:
		psi = [0] * len(x)
		for i in range(len(x)):
			for j,k,l in zip(q,p,C):
				psi[i] += l * ( norm )**0.25 * np.exp( -0.5 * m * omega * ( x[i] - j )**2 + 1j * k * ( x[i] - j ) )
	
		current_plot1.set_ydata( [ abs(_)**2 + e.real for _ in psi ] )
		txt = ax.text(np.sqrt(2 * e.real / m / omega**2) - 1,e.real+0.1,str(res.t+dt))
		plt.draw()
		plt.pause(1.5)
		ax.texts.remove(txt)
	else:
		print(res.t+dt,q[0].real,p[0].real,s.real,t.real,v.real,e.real,l.real)

