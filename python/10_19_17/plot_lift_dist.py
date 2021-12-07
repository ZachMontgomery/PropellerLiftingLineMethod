import time
import numpy as np
from propeller_mod import *
import matplotlib.pyplot as plt

start = time.time()
#prop1 = propeller('propeller.json')
prop2 = propeller('propeller.json')

x = np.zeros(prop2.nb)
x = prop2.cp[:prop2.nb,0]
tl = np.zeros(prop2.nb)
tg = np.zeros(prop2.nb)
dl = np.zeros(prop2.nb)
dg = np.zeros(prop2.nb)

tol_pitch = 1.e-5
tol_nm = 1.e-12
maxit = 100

J = float(input('Enter value for advance ratio: '))

#~ prop1.update_advance_ratio(J=J)
#~ # perform local iterative pitch solver
#~ prop1.local_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
#~ tl = -prop1.dF[:prop1.nb,2]
#~ dl = prop1.dF[:prop1.nb,1]

prop2.update_advance_ratio(J=J)
# perform global iterative pitch solver
prop2.global_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
tg = -prop2.dF[:prop2.nb,2]
dg = prop2.dF[:prop2.nb,1]
sg = prop2.dF[:prop2.nb,0]

f = open('pll_blade_loading_j{:3.1f}.txt'.format(J),'w')
for i in range(prop2.nb):
	line = '{0}		{1}		{2}		{3}\n'.format(prop2.cp[i,0], prop2.dF[i,1], prop2.dF[i,2], prop2.twist[i])
	f.write(line)
f.write('{}		{}'.format(prop2.CT, prop2.CL))

f.close()
close = time.time()
print('\n\nProgram ran for {:8.3f} seconds'.format(close - start))

plt.figure()
#~ plt.plot(x,tl)
plt.plot(x,tg)
plt.ylabel('Thrust (lbf)')
plt.xlabel('Spanwise length on blades (ft)')
#~ plt.legend(['Local Pitch Solver','Global Pitch Solver'])
plt.figure()
#~ plt.plot(x,dl)
plt.plot(x,dg)
plt.ylabel('Drag (lbf)')
plt.xlabel('Spanwise length on blades (ft)')
#~ plt.legend(['Local Pitch Solver','Global Pitch Solver'])
plt.figure()
#~ plt.plot(x,dl)
plt.plot(x,sg)
plt.ylabel('Side (lbf)')
plt.xlabel('Spanwise length on blades (ft)')
#~ plt.legend(['Local Pitch Solver','Global Pitch Solver'])


plt.show()
