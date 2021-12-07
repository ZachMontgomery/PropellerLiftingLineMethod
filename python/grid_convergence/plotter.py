import numpy as np
import re
# import json

numeric_const_pattern = r"""
     [-+]? # optional sign
     (?:
         (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
         |
         (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
     )
     # followed by optional exponent part if desired
     (?: [Ee] [+-]? \d+ ) ?
     """
rx = re.compile(numeric_const_pattern, re.VERBOSE)
del numeric_const_pattern

def myfunc(f):
	dummy = rx.findall(f.readline())
	return dummy

f = open('C_vs_J.txt', 'r')

zeta = np.zeros(11)
CT = np.zeros((11,7))
CP = np.zeros((11,7))
CL = np.zeros((11,7))
nodes = np.zeros(7)


for i in range(7):
	for j in range(11):
		
		zeta[j], CT[j,i], CP[j,i], CL[j,i], nodes[i] = myfunc(f)
		

import matplotlib.pyplot as plt


plt.figure()
plt.semilogx(nodes,CT[0,:],'o')
plt.title('J = 1.0')
plt.figure()
plt.semilogx(nodes,CT[1,:],'o')
plt.title('J = 0.9')
plt.figure()
plt.semilogx(nodes,CT[2,:],'o')
plt.title('J = 0.8')
plt.figure()
plt.semilogx(nodes,CT[3,:],'o')
plt.title('J = 0.7')
plt.figure()
plt.semilogx(nodes,CT[4,:],'o')
plt.title('J = 0.6')
plt.figure()
plt.semilogx(nodes,CT[5,:],'o')
plt.title('J = 0.5')
plt.figure()
plt.semilogx(nodes,CT[6,:],'o')
plt.title('J = 0.4')
plt.figure()
plt.semilogx(nodes,CT[7,:],'o')
plt.title('J = 0.3')
plt.figure()
plt.semilogx(nodes,CT[8,:],'o')
plt.title('J = 0.2')
plt.figure()
plt.semilogx(nodes,CT[9,:],'o')
plt.title('J = 0.1')
plt.figure()
plt.semilogx(nodes,CT[10,:],'o')
plt.title('J = 0.0')

plt.figure()
plt.plot(zeta,CT[:,1],'o')
plt.plot(zeta,CT[:,2],'o')
plt.plot(zeta,CT[:,3],'o')
plt.plot(zeta,CT[:,4],'o')
plt.plot(zeta,CT[:,5],'o')
plt.plot(zeta,CT[:,6],'o')


plt.show()
