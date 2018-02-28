import numpy as np
import matplotlib.pyplot as plt

def _norm_distrib(x,mu,sigma):
	return np.exp( -0.5 * ( (x - mu) / sigma )**2 ) / np.sqrt(2 * np.pi * sigma**2)

def _integration(x, y):
	integral = 0.0
	for i in range(len(x)):
		if i > 0:
			integral += 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1])
	return integral


N = 100000
x = [ 10 * np.random.random() - 5 for i in range(N)]
y = [ np.random.normal(0,1) for i in range(N)]
#plt.plot(x,[_norm_distrib(_,0,1) for _ in x])
#plt.show()
integral1 = sum( [ 3 * _**2 * _norm_distrib(_,0,1)  for _ in x ] ) * 10.0 / N
integral2 = sum( [ 3 * _**2  for _ in y ] ) * 1.0 / N
integral3 = _integration(np.linspace(-5,5,100000),[3 * _**2 * _norm_distrib(_,0,1) for _ in np.linspace(-5,5,100000)])
print(integral1, integral2, integral3)
