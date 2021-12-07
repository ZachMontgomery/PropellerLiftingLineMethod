import numpy as np
import matplotlib.pyplot as plt
import helix_func as hf

res = 100

x = np.zeros(res)
y = np.zeros(res)
yl = np. zeros(res)

yc = 1.
ymax = 10.

p = 1.5

plt.figure()


for i in range(res):
	
	x[i] = hf.interpolate(0., float(res-1), float(i), -1., 1.)
	
	
	if x[i] <= 0.:
		y[i] = yc * (1. - abs ( (-x[i]) ** p ) )
	else:
		y[i] = (ymax - yc) * x[i] ** (p+2.) + yc
	
	yl = y[i] * np.ones(2)
	
	plt.plot(np.array([-1.,1.]),yl,'r')


plt.plot(x, y, 'bs')
plt.show()

print(y)
