import numpy as np
import blade_shape as bs
import matplotlib.pyplot as plt

res = 1000

cl = np.zeros(res)
cla = np.zeros(res)
aoa = np.zeros(res)

for i in range(res):
	aoa[i] = 90. - float(i) / float(res-1) * 100.
	cl[i], cla[i] = bs.Cl(aoa[i] * np.pi / 180.)

plt.figure()
plt.plot(aoa,cl)
plt.plot(aoa,cla)
plt.show()
