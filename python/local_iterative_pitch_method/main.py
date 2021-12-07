import time
import numpy as np
from propeller_mod import *
import matplotlib.pyplot as plt

start = time.time()
prop1 = propeller('input.txt')
prop2 = propeller('input.txt')

res = 11
J = np.zeros(res)
CT1 = np.zeros(res)
CP1 = np.zeros(res)
CT2 = np.zeros(res)
CP2 = np.zeros(res)
tol_pitch = 1.e-8
tol_nm = 1.e-12
maxit = 1000
plt.ion()

for i in range(res):
	# renew advance ratio and update setup variables
	J[i] = 1. - float(i) / float(res-1) * 1.
	Ct_guess = .104 - .104 * J[i] ** 2. #float(input('Advance Ratio is {}, enter guess for CT: '.format(J[i])))
	
	prop1.update_advance_ratio(J=J[i],Ct=Ct_guess)
	
	# perform local iterative pitch solver
	prop1.local_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
	CT1[i] = prop1.CT
	CP1[i] = prop1.CP
	
	prop2.update_advance_ratio(J=J[i],Ct=Ct_guess)
	# perform global iterative pitch solver
	prop2.global_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
	CT2[i] = prop2.CT
	CP2[i] = prop2.CP
	
	k = i + 1
	if k == 1: k = res
	p1 = plt.figure(1)
	plt.clf()
	plt.plot(J[:k],CT1[:k])
	plt.plot(J[:k],CT2[:k])
	plt.ylabel('Coefficient of Thrust')
	plt.xlabel('Advance Ratio')
	plt.legend(['Local Pitch Solver','Global Pitch Solver'])
	p2 = plt.figure(2)
	plt.clf()
	plt.plot(J[:k],CP1[:k])
	plt.plot(J[:k],CP2[:k])
	plt.ylabel('Coefficient of Power')
	plt.xlabel('Advance Ratio')
	plt.legend(['Local Pitch Solver','Global Pitch Solver'])
	plt.draw()
	plt.pause(.1)

close = time.time()
print('\n\nProgram ran for {:8.3f} seconds'.format(close - start))

plt.figure(1)
plt.clf()
plt.plot(J,CT1)
plt.plot(J,CT2)
plt.ylabel('Coefficient of Thrust')
plt.xlabel('Advance Ratio')
plt.legend(['Local Pitch Solver','Global Pitch Solver'])
plt.figure(2)
plt.clf()
plt.plot(J,CP1)
plt.plot(J,CP2)
plt.ylabel('Coefficient of Power')
plt.xlabel('Advance Ratio')
plt.legend(['Local Pitch Solver','Global Pitch Solver'])
plt.show()
input('Program finished! Press [enter] to end.')
