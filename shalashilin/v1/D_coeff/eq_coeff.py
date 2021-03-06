import numpy as np
import matplotlib.pyplot as plt

def lagrangian(q,p,m,omega):
	return 0.5 * ( p**2 / m - m * omega**2 * q**2 )

def energy(q,p,m,omega):
	return 0.5 * ( p**2 / m + m * omega**2 * q**2 )

def func(t,q,p,m,omega):
	return 0.5 * ( omega * t + lagrangian(q,p,m,omega) * np.sin( 2 * omega * t ) / omega - 2 * q * p * np.sin( omega * t )**2 )  

def exponent(t,q,p,c,m,omega):
	return np.log(c) - 1j * 0.5 * ( omega * t + lagrangian(q,p,m,omega) * np.sin( 2 * omega * t ) / omega - 2 * q * p * np.sin( omega * t )**2 )

m = 1.0
omega = 1.0
A = 1.1584

#q = [A, 0.0, -A, 0.0]
#p = [0.0, A, 0.0, -A]

q = [A, -A, -A, A]
p = [A, A, -A, -A]

c = [1.0,1.0,1.0,1.0]

if len(q) == len(p) and len(p) == len(c):
	N = len(q)
else:
	raise NameError('inconsistent sizes of q, p, c')

t = np.linspace(0,20,1000)
fig = plt.figure( figsize = (10,6), facecolor = 'w' )
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
#fig1 = plt.figure( figsize = (10,6), facecolor = 'w' )
#ax3 = fig1.add_axes([0.1,0.1,0.8,0.8])
for i in range(N):
	e = [np.exp( exponent(_,q[i],p[i],c[i],m,omega) ) for _ in t]
	f = [func(_,q[i],p[i],m,omega) for _ in t]
	ax1.plot(t,[_.real for _ in e],label='q={:.4f}, p={:.4f}, C={:.4f}'.format(q[i],p[i],c[i]))
	ax1.set_title('Real part')
	ax2.plot(t,[_.imag for _ in e],label='q={:.4f}, p={:.4f}, C={:.4f}'.format(q[i],p[i],c[i]))
	ax2.set_title('Imaginary part')
#	ax3.plot(t,[_ for _ in f])

ax1.legend(loc="upper right")
ax2.legend(loc="upper right")

ax1.plot(t,[0 for _ in t])
ax2.plot(t,[0 for _ in t])
plt.show()


