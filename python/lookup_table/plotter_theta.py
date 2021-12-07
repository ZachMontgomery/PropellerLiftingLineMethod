import numpy as np
import matplotlib.pyplot as plt
import helix_func as hf

res = 360

x = np.zeros(res)
y = np.zeros(res)
yl = np. zeros(2)

ymax = 360.

p = 1.5

plt.figure()


for i in range(res):
	
	x[i] = float(i)
	
	y[i] = 360. * (1. - np.cos(float(i) * np.pi / float(res-1)))
	
	yl = y[i] * np.ones(2)
	
	plt.plot(np.array([0.,float(res)]),yl,'r')


plt.plot(x, y, 'bs')
plt.show()

print(y)
