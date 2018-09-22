import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
import fort_subroutines as fs

def _wave_packet( x, NumG, xi, eta, omega, D ):
	psi = [0]*len(x)
	for _ in range(len(x)):
		for i in range(NumG):
#			xi = m * omega[i] * q[i] + 1j * p[i]
#			eta = 0.25 * np.log( m * omega[i] / np.pi ) - 0.5 * ( m * omega[i] * q[i]**2 + 1j * q[i] * p[i] )
##			eta = 0.25 * np.log( m * omega[i] / np.pi ) - 0.5 * m * omega[i] * q[i]**2
			psi[_] += D[i] * np.exp(-0.5 * m * omega[i] * x[_]**2 + xi[i] * x[_] + eta[i])
#
	return psi


fin_ev =  open('../eigen_vectors.out','r')
fin_qp = open('../qp_init.out','r')

junk = [ fin_qp.readline() for i in range(2) ]
m = np.float_( fin_qp.readline() )

NumG = int(fin_ev.readline())

energy = np.float_( fin_ev.readline().split() )

with open('../init_cond.txt','r') as fin:
	junk = [ fin.readline() for i in range(2)]
	omega = [ float( fin.readline().split()[0].replace('D','E') ) ] * NumG

phase = [ 0.25 * np.log( x / np.pi ) for x in omega ]

vector = [ 0 for i in range(NumG)]
D = [ [0] for i in range(NumG) ] 
q = [ 0 for i in range(NumG) ]
p = [ 0 for i in range(NumG) ]  

for i in range(NumG):
	vector[i] = np.float_(fin_ev.readline().split() ).reshape(NumG,2)
	D[i] = [ vector[i][_][0] + 1j * vector[i][_][1] for _ in range(NumG) ]
	q[i], p[i], junk = np.float_( fin_qp.readline().split() )

fin_ev.close()
fin_qp.close()

D = np.transpose(D)
x = np.linspace(-5.5,5.5,500)
psi = [0]*NumG

fig_states = plt.figure(figsize = (6,6), facecolor = 'w' )
ax_states = fig_states.add_axes([0.1,0.1,0.85,0.85])
for i in range(NumG):
#	pprint(D[i])
	xi, eta = fs.change_var( q, p, omega, phase )
	psi[i] = _wave_packet(x, NumG, xi, eta, omega, D[i])
#	psi[i] = _wave_packet( x, NumG, m, omega, q, p, D[i] )
	ax_states.plot(x, [ abs(psi[i][_])**2 + energy[i] for _ in range(len(x)) ])
	ax_states.text( 4, i, energy[i] )
#	sumD = 0.0
#	for _ in range(NumG):
#		sumD += D[i][_]**2
#	D[i] /= np.sqrt(sumD)


#fig_weights = plt.figure(figsize = (6,6), facecolor = 'w')
#ax_w = fig_weights.add_axes([0.1,0.1,0.85,0.85])
#legend = [0] * NumG
#text = 'N: q, p\n'

#	ax_w.scatter( [ _+1 for _ in range(NumG) ], #
#    		      [ energy[i] for _ in range(NumG) ], # 
#		      s = [ 500 * abs(D[i][_])**2 for _ in range(NumG) ] )

#ax_w.scatter( q, [p[i]**2 * 0.5 + q[i]**2 * 9.072045 * 10**-7 for i in range(NumG)], s = [ 10 * abs(D[1][_])**2 for _ in range(NumG) ])

#	text += '{0}: {1}, {2}\n'.format(i+1,q[i],p[i])

#ax_states.plot(x, [ 9.072045 * 10**-7 * _**2 for _ in x])
ax_states.plot(x, [ 0.5 * _**2 for _ in x])

#x_coord = NumG + 1
#ax_w.text(x_coord, 0.005, text, fontsize = 15)
#ax_w.set_ylim([-0.0001,24.])
plt.show()


