import numpy as np
import matplotlib.pyplot as plt

def lagrangian(q,p,m,omega):
	return 0.5 * ( p**2 / m - m * omega**2 * q**2 )

def exponent(t,q,p,c,m,omega):
	return np.log(c) - 1j * 0.5 * ( omega * t + lagrangian(q,p,m,omega) * np.sin( 2 * omega * t ) / omega)

m = 1.0
omega = 1.0
A = 2.09336
#A = 1.8315
q = [-A, 0.0, 0.0, A]
p = [0.0, -A, A, 0.0]
c = [1.0,1.0,1.0,1.0]

if len(q) == len(p) and len(p) == len(c):
	N = len(q)
else:
	raise NameError('inconsistent sizes of q, p, c')

t = np.linspace(0,20,1000)
fig = plt.figure( figsize = (10,6), facecolor = 'w' )
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
for i in range(N):
	e = [np.exp( exponent(_,q[i],p[i],c[i],m,omega) ) for _ in t]
	ax1.plot(t,[_.real for _ in e],label='q={:.4f}, p={:.4f}, C={:.4f}'.format(q[i],p[i],c[i]))
	ax1.set_title('Real part')
	ax2.plot(t,[_.imag for _ in e],label='q={:.4f}, p={:.4f}, C={:.4f}'.format(q[i],p[i],c[i]))
	ax2.set_title('Imaginary part')

ax1.legend(loc="upper right")
ax2.legend(loc="upper right")

ax1.plot(t,[0 for _ in t])
ax2.plot(t,[0 for _ in t])
plt.show()


