from propeller_mod import *
import numpy as np


prop = propeller('input.txt')


res = 50
J = np.zeros(res)
CL = np.zeros(res)
CD = np.zeros(res)
f = open('comparison.txt', 'w')

CL_machup = np.ones(res) * .8665
CD_machup = np.ones(res) * .0286

for j in range(res):
	if j <= res / 2:
		J[j] = 10. + float(j) / float(res/2) * 90.
	else:
		J[j] = float(j-res/2) / float(res/2-1) * 900. + 100.
	#J[j] = 50.
	#prop.vinf = J[j] * prop.w * prop.rp / np.pi
	prop.w = prop.vinf * np.pi / prop.rp / J[j]
	#prop.Vinf[2] = prop.vinf
	prop.b = prop.vinf / prop.w * 2. * np.pi
	prop.W[2] = prop.rl * prop.w
	
	print('Freestream is {} and pitch is {} and omega is {}'.format(prop.vinf,prop.b,prop.w))
	
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
		
		print('Advance ratio is {}'.format(J[j]))
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
	
	CL[j] = prop.total_force[1] * 2. / prop.rho / prop.vinf ** 2. / prop.c[0,0] / (prop.rp - prop.rh)
	CD[j] = prop.total_force[2] * 2. / prop.rho / prop.vinf ** 2. / prop.c[0,0] / (prop.rp - prop.rh)
	
	line = '{0}		{1}		{2}\n'.format(J[j], CL[j], CD[j])
	f.write(line)
	print(line)

f.close()

import matplotlib.pyplot as plt

plt.figure()
plt.plot(J, CL)
plt.plot(J, CL_machup)

plt.figure()
plt.plot(J, CD)
plt.plot(J, CD_machup)

#plt.figure()
#plt.plot(J, CP)

plt.show()


