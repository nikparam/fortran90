import os
import matplotlib.pyplot as plt
import numpy as np
import time
import potential

def _wave_packet( x, NumG, m, omega, q, p, D ):
	psi = [0]*len(x)
	for _ in range(len(x)):
		for i in range(NumG):
			xi = m * omega[i] * q[i] + 1j * p[i]
#			eta = 0.25 * np.log( 0.5 * m * omega[i] / np.pi ) - 0.5 * ( m * omega[i] * q[i]**2 + 1j * q[i] * p[i] )
			eta = 0.25 * np.log(  m * omega[i] / np.pi ) - 0.5 * m * omega[i] * q[i]**2
			psi[_] += D[i] * np.exp(-0.5 * m * omega[i] * x[_]**2 + xi * x[_] + eta)

	return psi

def _ffmt( i ):
	if i < 10:
		return '0' + str(i)
	elif i < 100:
		return str(i)

def _read( fname ):
	return np.float_( [string.replace('D','E') for string in fname.readline().split()] )

def _morse( x ):
	return 0.05 * (1 - np.exp( -0.005 * (x - 200) ))**2


BASEDIR = os.getcwd()
path_all = BASEDIR + '/all/'
path_potential = BASEDIR + '/potential/'
init_cond_file = open(path_all + 'init_cond.txt','r')

switch = np.float(init_cond_file.readline().split()[0])

NumG, Tmax, Tstep = _read( init_cond_file )

m = _read( init_cond_file )
dE = _read( init_cond_file )

omega = [np.float(init_cond_file.readline().split()[0].replace('D','E')) for _ in range(int(NumG))]

x = np.linspace(-3.5,3.5,200)

start_time = time.time()

NumSteps = int(Tmax / Tstep)

coeff_fin = open(path_all+'D_coeff.out', 'r')
qp_fin = open(path_all + 'coords.out', 'r')

plt.ion()
ax = plt.gca()
ax.set_autoscale_on(True)
ax.set_ylim([1,3])
psi = [0]*len(x)
ax.plot(x, [ _**2 for _ in x])

current_plot1, = ax.plot( x, [abs(psi[_])**2 for _ in range(len(psi))] )
count = 0
num = 0
sum_psi = [ 0 ] * len(x)
NumG = int(NumG)
for line in coeff_fin:
	D_ufmt = [0] * (2 * NumG)
	D_norm = [0] * NumG
	q = [0] * NumG
	p = [0] * NumG
	try:
		t = np.float( line.split()[0] )
		D_ufmt = np.float_( line.split()[1:2*NumG+1] ).reshape((NumG,2))
		D_fmt = [ D_ufmt[i][0] + 1j * D_ufmt[i][1] for i in range(NumG) ]
		D_norm = np.float_( line.split()[2*NumG+1:3*NumG+1] )
		D_total_norm = np.float( line.split()[3*NumG+1] )
		qp_line = qp_fin.readline()

		for i in range(NumG):
			q[i] = np.float( qp_line.split()[i+1] )
			p[i] = np.float( qp_line.split()[-1] )

	except:
		pass


	if bin(count)[-13:] == '0'*13 or count == 0:
		psi = _wave_packet(x, NumG, m, omega, q, p, D_fmt)
		current_plot1.set_ydata([ abs(psi[i])**2 + dE + 0.5  for i in range(len(psi))])
		ax.autoscale_view(True,True)
		txt = ax.text(3.0, 0.075, str(t))
		plt.draw()
		plt.pause(1)
		ax.texts.remove(txt)
		num += 1
	count += 1

plt.show()
print('CPU time= {:f} '.format( time.time() - start_time ) )
