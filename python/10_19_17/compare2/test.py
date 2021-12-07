import time

from propeller_mod import *
import numpy as np

start = time.time()



prop = propeller('propeller.json')


res = 1
J = np.zeros(res)
CT = np.zeros(res)
CL = np.zeros(res)
CP = np.zeros(res)
#f = open('pll_C_vs_J.txt', 'w')

for j in range(res):
	#J[j] = float(j) / float(res-1)
	J[j] = 0.4
	
	
	prop.global_iterative_pitch_solver(1.e-5,1.e-12,200)
	
	J[j] = prop.J
	CT[j] = prop.CT
	CL[j] = prop.CL
	CP[j] = prop.CP
	
	line = '{0}		{1}		{2}		{3}\n'.format(J[j], CT[j], CP[j], CL[j])
	#f.write(line)

#f.close()
close = time.time()

print('Program ran for {:8.3f} seconds'.format(close - start))

#import matplotlib.pyplot as plt

#plt.figure()
#plt.plot(J, CT)

#plt.figure()
#plt.plot(J, CL)

#plt.figure()
#plt.plot(J, CP)

#plt.show()


