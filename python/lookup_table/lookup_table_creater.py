import numpy as np
import helix_func as hf

res_theta = 360
res_sigma = 100
res_lambda = 700

lookup_table = np.zeros((res_lambda,res_sigma,res_theta,3))
count = 0
cyl = np.zeros(3)
blah = np.zeros(3)

theta_low = 0.
theta_high = 2. * np.pi

sigma_low = 0.
sigma_high = 10.

lambda_low = .001
lambda_high = 70.

n_low = 500.
n_high = 20.

m_low = 20.
m_high = 500.

p_low = 1.
p_high = 3.

rl = -1.

for j in range(res_lambda):
	
	lamb = hf.interpolate(0.,float(res_lambda-1),float(j),lambda_low,lambda_high)
	#print('	lambda is ',lamb)
	
	# determine n, m, p
	n = hf.interpolate(0.,float(res_lambda-1),float(j),n_low,n_high)
	#print('	n is ',n)
	m = hf.interpolate(0.,float(res_lambda-1),float(j),m_low,m_high)
	#print('	m is ',m)
	p = hf.interpolate(0.,float(res_lambda-1),float(j),p_low,p_high)
	#print('	p is ',p)
	#print('\n\nend of helix is {}'.format(lamb*n))
	#input()
	
	for k in range(res_sigma):
		
		zeta = hf.interpolate(0., float(res_sigma-1), float(k), -1., 1.)
		if zeta <= 0.:
			sigma = 1. - (-zeta) ** 1.5
		else:
			sigma = (sigma_high - 1.) * zeta ** 3.5 + 1.
		cyl[0] = sigma
		print('		sigma is {} and k is {}'.format(cyl[0],k))
		
		
		for l in range(res_theta):
			
			theta = theta_high * (1. - np.cos(float(l) * np.pi / float(res_theta-1)))
			#print('			theta is ',theta)
			
			if cyl[0] == 1. and (theta == 0. or theta == 2.*np.pi):
				print('Point lies on helix')
				print('coordinates: {}:{}'.format(k,l))
				input()
				
			
			cyl[1] = theta
			
			
			#lookup_table[j,k,l,:] = hf.integrate_nondim(blah, lamb, rl, cyl, n, m, p, 0., 'r')
			count += 1
			
		
		#print('{:03.5f}% done'.format(  float(count)/float(res_theta*res_sigma*res_lambda)*100.  ))
	



#np.save('lookup_table',lookup_table)

#print(lookup_table)
