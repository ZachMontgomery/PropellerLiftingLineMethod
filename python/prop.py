import time
import numpy as np
from propeller_mod import *

start = time.time()

prop = propeller('input.txt')

#J[j] = 0.4
#prop.vinf = J[j] * prop.w * prop.rp / np.pi
#prop.Vinf[2] = prop.vinf
#prop.b = prop.vinf / prop.w * 2. * np.pi

#prop.renew_Vp()

print('Advance ratio is {}\n\n'.format(prop.J))



count = 1
ea_max = 1.

while ea_max > 1.e-5:
	
	print('\n\n\nPitch Iteration {}'.format(count))
	#print('The helical pitch is {}\n\n\n'.format(prop.b))
	
	#print(prop.b)
	#input()
	
	prop.renew_nu()
	prop.pll(1.e-12, 100)
	
	# individual pitch solvers per node
	temp = np.zeros([prop.cp_total,2])
	for i in range(prop.nb):
		for j in range(prop.n_blades):
			k = i + j * prop.nb
			temp[k,0] = prop.V[k,2] / prop.w * 2. * np.pi
			if i < prop.nb-1:
				temp[k,1] = prop.V[k+1,2] / prop.w * 2. * np.pi
			else:
				temp[k,1] = temp[k,0]
	
	ea_max = 0.
	for i in range(prop.cp_total):
		for j in range(2):
			ea = abs((temp[i,j] - prop.b[i,j]) / temp[i,j])
			if ea > ea_max: ea_max = ea
	print('Max Approximate Error for the pitch solver is {}%\n\n'.format(ea_max*100.))
	count += 1
	
	for i in range(prop.cp_total):
		for j in range(2):
			prop.b[i,j] += prop.Omega1 * (temp[i,j] - prop.b[i,j])
	
	
	
	# calculate average axial velocity across propeller blades
	#avg_v = 0.
	#for i in range(prop.cp_total):
	#	avg_v += prop.V[i,2]
	#avg_v = avg_v / float(prop.cp_total)
	# calculate how long the helix pitch should be based on the average axial velocity
	#temp = avg_v / prop.w * 2. * np.pi
	
	#print('The helical pitch length due to the induced velocity should be {}'.format(temp))
	#ea = abs((temp - prop.b) / temp)
	#print('Approximate error for the pitch solver is {}%\n\n'.format(ea*100.))
	#count += 1
	# update pitch length
	#prop.b += prop.Omega1 * (temp - prop.b)
	
	


close = time.time()
print('\n\nProgram ran for {:8.3f} seconds'.format(close - start))
