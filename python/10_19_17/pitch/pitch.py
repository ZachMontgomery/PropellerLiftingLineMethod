import numpy as np
import matplotlib.pyplot as plt

res1 = 100
res2 = 10

zeta = np.linspace(0.1, 1.0, num=res1)
beta = np.zeros((res1,res2))

tanal0 = np.tan(-2.1 * np.pi / 180.)

for i in range(res2):
	
	kc = float(i) / float(res2 - 1) * 1.
	
	beta[:,i] = np.arctan( (kc - np.pi * zeta * tanal0) / (np.pi * zeta + kc * tanal0) ) * 180. / np.pi

plt.plot(zeta, beta)
plt.show()
