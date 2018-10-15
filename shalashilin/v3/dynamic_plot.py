import os
import matplotlib.pyplot as plt
import numpy as np
import time
import potential

def _wave_packet( x, NumG, m, omega, q, p, D ):
	psi = [0]*len(x)
	for _ in range(len(x)):
		for i in range(NumG):
			psi[_] += D[i] * np.exp(-0.5 * m * omega[i] * ( x[_] - q[i] )**2 + 1j * p[i] *  x[_] )
	return psi

def _ffmt( i ):
	if i < 10:
		return '0' + str(i)
	elif i < 100:
		return str(i)

def _read( fname ):
	return [ string.replace('D','E') for string in fname.readline().split() ]

def _morse( x ):
	return 0.05 * (1 - np.exp( -0.005 * (x - 200) ))**2


BASEDIR = os.getcwd()
path_all = BASEDIR + '/all/'
path_potential = BASEDIR + '/potential/'
init_cond_file = open(path_all + 'init_cond.txt','r')

mean = init_cond_file.readline().split()[0]
switch = np.float(init_cond_file.readline().split()[0])
lambda_r = _read( init_cond_file )[0]

NumG, Tmax, Tstep, = _read( init_cond_file )[0:3]
dE = _read( init_cond_file )[0]
m = _read( init_cond_file )[0]

NumG = int(NumG)
Tmax = float(Tmax)
Tstep = float(Tstep)
dE = float(dE)
m = float(m)

omega = [np.float(init_cond_file.readline().split()[0].replace('D','E')) for _ in range(int(NumG))]

x = np.linspace(-5,40,200)

start_time = time.time()

NumSteps = int(Tmax / Tstep)

coeff_fin = open(path_all+'D_coeff.out', 'r')
qp_fin = open(path_all + 'coords.out', 'r')

plt.ion()
ax = plt.gca()
ax.set_autoscale_on(True)
ax.set_ylim([0,5])
psi = [0]*len(x)
#ax.plot(x, [ 0.25 *  _**2 for _ in x])
ax.plot(x, [np.exp(-_) for _ in x])

#current_plot1, current_plot2, current_plot3, current_plot4 = ax.plot( x, [abs(psi[_])**2 for _ in range(len(psi))], #
#								      x, [abs(psi[_])**2 for _ in range(len(psi))], #
#								      x, [abs(psi[_])**2 for _ in range(len(psi))], #
#								      x, [abs(psi[_])**2 for _ in range(len(psi))] )

current_plot4, = ax.plot( x, [abs(psi[_])**2 for _ in range(len(psi))])

count = 0
num = 0
sum_psi = [ 0 ] * len(x)
NumG = int(NumG)
for line in coeff_fin:
	D_ufmt = [0] * (2 * NumG)
	D_norm = [0] * NumG
	q = [0] * NumG
	p = [0] * NumG

	t = np.float( line.split()[0] )
	D_ufmt = np.float_( line.split()[1:2*NumG+1] ).reshape((NumG,2))
	D_fmt = [ D_ufmt[i][0] + 1j * D_ufmt[i][1] for i in range(NumG) ]
	D_norm = np.float_( line.split()[2*NumG+1:3*NumG+1] )
	qp_line = qp_fin.readline()

	for i in range(NumG):
		q[i] = np.float( qp_line.split()[i+1] )
		p[i] = np.float( qp_line.split()[i+NumG+1] )

#	psi1 = _wave_packet(x, 1, m, omega, [ q[0] ], [ p[0] ], [ D_fmt[0] ])
#	psi2 = _wave_packet(x, 1, m, omega, [ q[1] ], [ p[1] ], [ D_fmt[1] ])
#	psi3 = _wave_packet(x, 1, m, omega, [ q[2] ], [ p[2] ], [ D_fmt[2] ])
	psi4 = _wave_packet(x, NumG, m, omega, q, p, D_fmt)
#	current_plot1.set_ydata([ abs( psi1[i] ) / abs(D_fmt[0])  for i in range(len(psi))])
#	current_plot2.set_ydata([ abs( psi2[i] ) / abs(D_fmt[1])  for i in range(len(psi))])
#	current_plot3.set_ydata([ abs( psi3[i] ) / abs(D_fmt[2])  for i in range(len(psi))])
	current_plot4.set_ydata([ abs( psi4[i] )  for i in range(len(psi))])
	ax.autoscale_view(True,True)
	txt = ax.text(3.0, 0.075, str(t))
	plt.draw()
	plt.pause(0.1)
	ax.texts.remove(txt)
	num += 1

plt.show()
print('CPU time= {:f} '.format( time.time() - start_time ) )
