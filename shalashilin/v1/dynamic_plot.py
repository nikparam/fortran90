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
			eta = 0.25 * np.log( 0.5 * m * omega[i] / np.pi ) - 0.5 * ( m * omega[i] * q[i]**2 + 1j * q[i] * p[i] )
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
path_traj = BASEDIR + '/mpi/'
path_coeff = BASEDIR + '/D_coeff/'
path_potential = BASEDIR + '/potential/'
init_cond_file1 = open(path_traj + 'init_cond.txt','r')
init_cond_file2 = open(path_coeff + 'init_cond.txt','r')

NumG = [0]*2
Tmax = [0]*2
Tstep = [0]*2
m = [0]*2
switch = [0]*2

switch[0] = np.float(init_cond_file1.readline())
switch[1] = np.float(init_cond_file2.readline())

NumG[0], Tmax[0], Tstep[0] = _read( init_cond_file1 )
NumG[1], Tmax[1], Tstep[1] = _read( init_cond_file2 )

m[0] = _read( init_cond_file1 )
dE = _read( init_cond_file2 )
m[1] = _read( init_cond_file2 )

omega = [np.float(init_cond_file2.readline().split()[0].replace('D','E')) for _ in range(int(NumG[1]))]

x = np.linspace(-2,20,200)
start_time = time.time()
if NumG[0] == NumG[1] and Tstep[0] == Tstep[1] and m[0] == m[1] and switch[0] == switch[1]:
	switch = int(switch[0])
	NumG = int(NumG[1])
	Tstep = float(Tstep[1])
	m = float(m[1])
	Tmax = float(min(Tmax))
	NumSteps = int(Tmax / Tstep)

	coeff_fin = open(path_coeff+'D_coeff.out','r')
	qp_fin = [ open(path_traj + 'out' + _ffmt(_) + '.out','r') for _ in range(1,NumG+1) ]
#	params = potential.read_params(switch,'./potential/params.txt')
#	print(params)

	plt.ion()
	ax = plt.gca()
	ax.set_autoscale_on(True)
	ax.set_ylim([0,1.5])
	psi = [0]*len(x)
#	ax.plot(x, [potential.potential_energy(_,params) for _ in x])
#	ax.plot(x, [ ( 1 - np.exp(-i) )**2 for i in x ] )
	ax.plot(x, [ np.exp(-i) for i in x ] )
#	current_plot1, current_plot2, = ax.plot( x, [psi[_].real for _ in range(len(psi))], \
#					 	  x, [psi[_].imag for _ in range(len(psi))])

	current_plot1, = ax.plot( x, [abs(psi[_])**2 for _ in range(len(psi))] )
	count = 0
	num = 0
	sum_psi = [ 0 ] * len(x)
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
			for i in range(NumG):
				qp_line = qp_fin[i].readline()
				q[i] = np.float( qp_line.split()[1] )
				p[i] = np.float( qp_line.split()[2] )
		except:
			pass


		if bin(count)[-13:] == '0'*13 or count == 0:
#			print(q,p)
			psi = _wave_packet(x, NumG, m, omega, q, p, D_fmt)
#			sum_psi = [ sum_psi[i] + abs( psi[i] )**2 for i in range( len(x) ) ]
			current_plot1.set_ydata([ abs(psi[i])**2  for i in range(len(psi))])
#			current_plot1.set_ydata([ psi[i].real for i in range(len(psi))])
#			ax.relim()
			ax.autoscale_view(True,True)
			txt = ax.text(3.0, 0.075, str(t))
			plt.draw()
			plt.pause(1)
			ax.texts.remove(txt)
			num += 1
		count += 1

#	ax.plot( x, [ sum_psi[i] / num + 2.5 for i in range( len(x) ) ])	

plt.show()
print('CPU time= {:f} '.format( time.time() - start_time ) )
