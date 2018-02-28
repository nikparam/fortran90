import matplotlib.pyplot as plt
import numpy as np
from functools import partial
from time import time

D = [ [1, 1, 1, 1], [-1, 1, -1j, 1j], [1, 1, -1, -1], [ 1j, -1j, 1, -1] ] 
A = 0.88622693

g = lambda i, x: abs(  D[i][0] * np.exp( -0.5 * x**2 - 1j * A * x) + \
	               D[i][1] * np.exp( -0.5 * x**2 + 1j * A * x) + \
	               D[i][2] * np.exp( -0.5 * ( x + A )**2 ) +  \
		       D[i][3] * np.exp( -0.5 * ( x - A )**2 ) ) **2
a = 5
NPoints = 500
x = np.linspace(-a,a,NPoints)

df = lambda f, x1, dx: ( f( x1 + dx ) - f( x1 ) ) / dx

start_time = time()

for i in range(4):
	h = partial(g,i)
	dh = partial(df,h)
	plt.plot( x,list(map(h,x)) ) 

	step = 2 * a / NPoints
	df0 = 0.0
	df1 = 0.0
	count = 0
	for j in x:
		df0 = dh(j,step)
		if df1 < 0 and df0 > 0: 
			count += 1
#			print( 'x = {0:g}, psi[{1}](x) = {2:g}'.format( j , i, h(j) ) )
		df1 = df0

	print( 'psi[{0}] zero count = {1}'.format( i, count ) )

print( 'CPU time = {0}'.format( time() - start_time ) )

plt.show()
