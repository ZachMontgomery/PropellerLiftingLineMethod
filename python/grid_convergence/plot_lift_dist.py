import time
import numpy as np
from propeller_mod import *
import matplotlib.pyplot as plt




start = time.time()
prop = propeller('propeller.json')

x = np.copy(prop.cp[:prop.nb,0])
t = np.zeros(prop.nb)
d = np.zeros(prop.nb)

tol_pitch = 1.e-5
tol_nm = 1.e-12
maxit = 100


method = input('Use the local pitch solver (y/n): ')

# prop.update_advance_ratio(J=prop.J)

flag = 0
while flag == 0:
	if method == 'y':
		# perform local iterative pitch solver
		prop.local_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
		flag = 1
	elif method == 'n':
		# perform global iterative pitch solver
		prop.global_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
		flag = 1
	else:
		method = input('Wrong input! Enter \'y\' or \'n\' ')

t = -prop.dF[:prop.nb,2]
d = prop.dF[:prop.nb,1]
s = prop.dF[:prop.nb,0]

f = open('pll_blade_loading.txt','w')
for i in range(prop.nb):
	line = '{0}		{1}		{2}		{3}\n'.format(prop.cp[i,0], prop.dF[i,1], prop.dF[i,2], prop.twist[i])
	f.write(line)
f.write('{}		{}'.format(prop.CT, prop.CL))

f.close()
close = time.time()
print('\n\nProgram ran for {:8.3f} seconds'.format(close - start))

plt.figure()
plt.plot(x,t)
plt.ylabel('Thrust (lbf)')
plt.xlabel('Spanwise length on blades (ft)')
plt.figure()
plt.plot(x,d)
plt.ylabel('Drag (lbf)')
plt.xlabel('Spanwise length on blades (ft)')
plt.figure()
plt.plot(x,s)
plt.ylabel('Side (lbf)')
plt.xlabel('Spanwise length on blades (ft)')
plt.show()

