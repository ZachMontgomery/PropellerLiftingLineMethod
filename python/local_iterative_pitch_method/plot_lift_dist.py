import time
import numpy as np
from propeller_mod import *
import matplotlib.pyplot as plt

start = time.time()
prop1 = propeller('input.txt')
prop2 = propeller('input.txt')

x = np.zeros(prop1.nb)
x = prop1.cp[:prop1.nb,0]
tl = np.zeros(prop1.nb)
tg = np.zeros(prop1.nb)
dl = np.zeros(prop1.nb)
dg = np.zeros(prop1.nb)

tol_pitch = 1.e-10
tol_nm = 1.e-12
maxit = 20

J = .25

prop1.update_advance_ratio(J=J)
# perform local iterative pitch solver
prop1.local_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
tl = -prop1.dF[:prop1.nb,2]
dl = prop1.dF[:prop1.nb,1]

prop2.update_advance_ratio(J=J)
# perform global iterative pitch solver
prop2.global_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
tg = -prop2.dF[:prop2.nb,2]
dg = prop2.dF[:prop2.nb,1]

close = time.time()
print('\n\nProgram ran for {:8.3f} seconds'.format(close - start))

plt.figure()
plt.plot(x,tl)
plt.plot(x,tg)
plt.ylabel('Thrust (lbf)')
plt.xlabel('Spanwise length on blades (ft)')
plt.legend(['Local Pitch Solver','Global Pitch Solver'])
plt.figure()
plt.plot(x,dl)
plt.plot(x,dg)
plt.ylabel('Drag (lbf)')
plt.xlabel('Spanwise length on blades (ft)')
plt.legend(['Local Pitch Solver','Global Pitch Solver'])
plt.show()
