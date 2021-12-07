from propeller_mod import *
import numpy as np


prop = propeller('input.txt')


res = prop.nb
zeta = np.zeros(res)
x = np.zeros(res)
y = np.zeros(res)
z = np.zeros(res)

J = .5
prop.vinf = J * prop.w * prop.rp / np.pi
prop.Vinf[2] = prop.vinf
prop.b = prop.vinf / prop.w * 2. * np.pi

prop.renew_Vp()


count = 0
ea = 1.

while ea > 1.e-10:
	
	prop.renew_nu()
	
	prop.pll(1.e-12, 100)
	
	avg_v = 0.
	
	for i in range(prop.n_blades * prop.nb):
		avg_v += prop.V[i,2]
	
	avg_v = avg_v / float(prop.n_blades * prop.nb)
	
	temp = avg_v / prop.w * 2. * np.pi
	
	print('Advance ratio is {}'.format(J))
	#print('The helical pitch length due to the induced velocity should be {}'.format(temp))
	#print('The helical pitch currently is {}'.format(prop.b))
	ea = abs((temp - prop.b) / temp)
	print('Approximate error for the pitch solver is {}%'.format(ea*100.))
	count += 1
	print('Iteration {}'.format(count))
	prop.b += prop.Omega1 * (temp - prop.b)
	
	print('The new helical pitch is {}'.format(prop.b))
	#input()
	
	
	#temp = prop.vinf / prop.w * 2. * np.pi
	
	#input()


f = open('pll_blade_loading.txt', 'w')

for i in range(prop.nb):
	zeta[i] = prop.cp[i,0]
	x[i] = prop.dF[i,0]
	y[i] = prop.dF[i,1]
	z[i] = prop.dF[i,2]
	
	line = '{0}		{1}		{2}		{3}\n'.format(prop.cp[i,0], prop.dF[i,0], prop.dF[i,1], prop.dF[i,2])
	f.write(line)

f.close()

import matplotlib.pyplot as plt

plt.figure()
plt.plot(zeta, x)
plt.plot(zeta, y)
plt.plot(zeta, -z)

plt.show()


