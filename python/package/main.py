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
ea = 1.

while ea > 1.e-10:
	
	print('\n\n\nPitch Iteration {}'.format(count))
	print('The helical pitch is {}\n\n\n'.format(prop.b))
	
	prop.renew_nu()
	prop.pll(1.e-12, 100)
	
	# calculate average axial velocity across propeller blades
	avg_v = 0.
	for i in range(prop.cp_total):
		avg_v += prop.V[i,2]
	avg_v = avg_v / float(prop.cp_total)
	# calculate how long the helix pitch should be based on the average axial velocity
	temp = avg_v / prop.w * 2. * np.pi
	
	print('The helical pitch length due to the induced velocity should be {}'.format(temp))
	ea = abs((temp - prop.b) / temp)
	print('Approximate error for the pitch solver is {}%\n\n'.format(ea*100.))
	count += 1
	# update pitch length
	prop.b += prop.Omega1 * (temp - prop.b)
	
	


close = time.time()
print('\n\nProgram ran for {:8.3f} seconds'.format(close - start))
