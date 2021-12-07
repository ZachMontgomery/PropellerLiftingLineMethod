import numpy as np
import helix_func as hf

lookup_table = np.load('lookup_table.npy')


def read_table_linear(a, b, x, y):
	lamb = b / a
	sigma = (x**2. + y**2.)**.5 / a
	theta = np.arctan2(y,x)
	
	v = np.zeros(3)
	
	theta_low = 0.
	theta_high = 2. * np.pi
	
	sigma_low = 0.
	sigma_high = 10.
	
	lambda_low = .001
	lambda_high = 70.
	
	res_theta = 360.
	res_sigma = 100.
	res_lambda = 700.
	
	#*******************************************************************
	#*******************************************************************
	#*******************************************************************
	#*******************************************************************
	if sigma <= 1.:
		zeta = -(1. - sigma) ** (2./3.)
	else:
		zeta = ((sigma - 1.) / (sigma_high - 1.)) ** (2./7.)
	k = hf.interpolate(-1.,1.,zeta,0.,float(res_sigma-1))
	
	
	l = float(res_theta-1) / np.pi * np.arccos(1. - theta / theta_high)
	#*******************************************************************
	#*******************************************************************
	#*******************************************************************
	#*******************************************************************
	
	
	j = hf.interpolate(lambda_low, lambda_high, lamb, 0., res_lambda-1.)
	
	j0 = int(np.floor(j))
	j1 = int(np.ceil(j))
	k0 = int(np.floor(k))
	k1 = int(np.ceil(k))
	l0 = int(np.floor(l))
	l1 = int(np.ceil(l))
	
	if j1 == j0: j1 += 1
	if k1 == k0: k1 += 1
	if l1 == l0: l1 += 1
	
	if k0 == 49 or k0 == 50 or k1 == 49 or k1 == 50:
		if sigma <= 1.:
			k1 = 49
			k0 = 48
		else:
			k0 = 50
			k1 = 51
	
	print(j0,j1,k0,k1,l0,l1,'indecies')
	
	jd = (j - float(j0)) / (float(j1) - float(j0))
	kd = (k - float(k0)) / (float(k1) - float(k0))
	ld = (l - float(l0)) / (float(l1) - float(l0))
	
	V000 = lookup_table[j0,k0,l0,:]
	V100 = lookup_table[j1,k0,l0,:]
	V110 = lookup_table[j1,k1,l0,:]
	V010 = lookup_table[j0,k1,l0,:]
	V001 = lookup_table[j0,k0,l1,:]
	V101 = lookup_table[j1,k0,l1,:]
	V111 = lookup_table[j1,k1,l1,:]
	V011 = lookup_table[j0,k1,l1,:]
	
	c00 = V000 * (1. - jd) + V100 * jd
	c01 = V001 * (1. - jd) + V101 * jd
	c10 = V010 * (1. - jd) + V110 * jd
	c11 = V011 * (1. - jd) + V111 * jd
	
	c0 = c00 * (1. - kd) + c10 * kd
	c1 = c01 * (1. - kd) + c11 * kd
	
	v = c0 * (1. - ld) + c1 * ld
	
	v = v / a
	
	return v

def read_table_quadratic(a1, b, x, y):
	lamb = b / a1
	sigma = (x**2. + y**2.)**.5 / a1
	theta = np.arctan2(y,x)
	
	v = np.zeros(3)
	
	theta_low = 0.
	theta_high = 2. * np.pi
	
	sigma_low = 0.
	sigma_high = 10.
	
	lambda_low = .001
	lambda_high = 70.
	
	res_theta = 360.
	res_sigma = 100.
	res_lambda = 700.
	
	#*******************************************************************
	#*******************************************************************
	#*******************************************************************
	#*******************************************************************
	if sigma <= 1.:
		zeta = -(1. - sigma) ** (2./3.)
	else:
		zeta = ((sigma - 1.) / (sigma_high - 1.)) ** (2./7.)
	k = hf.interpolate(-1.,1.,zeta,0.,float(res_sigma-1))
	
	
	l = float(res_theta-1) / np.pi * np.arccos(1. - theta / theta_high)
	#*******************************************************************
	#*******************************************************************
	#*******************************************************************
	#*******************************************************************
	
	
	j = hf.interpolate(lambda_low, lambda_high, lamb, 0., res_lambda-1.)
	
	j0 = int(np.floor(j))
	j1 = int(np.ceil(j))
	k0 = int(np.floor(k))
	k1 = int(np.ceil(k))
	l0 = int(np.floor(l))
	l1 = int(np.ceil(l))
	
	if j1 == j0: j1 += 1
	if k1 == k0: k1 += 1
	if l1 == l0: l1 += 1
	
	
	jd = (j - float(j0)) / (float(j1) - float(j0))
	kd = (k - float(k0)) / (float(k1) - float(k0))
	ld = (l - float(l0)) / (float(l1) - float(l0))
	
	if jd < .5:
		j2 = j0 - 1
		jp = 'l'
	else:
		j2 = j1 + 1
		jp = 'u'
	
	if k0 == 49 or k0 == 50 or k1 == 49 or k1 == 50:
		if sigma <= 1.:
			k2 = 49
			k1 = 48
			k0 = 47
			kp = 'u'
		else:
			k0 = 50
			k1 = 51
			k2 = 52
			kp = 'u'
	elif kd < .5:
		k2 = k0 - 1
		kp = 'l'
	else:
		k2 = k1 + 1
		kp = 'u'
	
	if l0 == 0:
		l2 = l1 + 1
		lp = 'u'
	elif l1 >= 359:
		l2 = l0 - 1
		lp = 'l'
	elif ld < .5:
		l2 = l0 - 1
		lp = 'l'
	else:
		l2 = l1 + 1
		lp = 'u'
	
	print(j0,j1,j2)
	print(k0,k1,k2)
	print(l0,l1,l2,'indecies')
	
	V000 = lookup_table[j0,k0,l0,:]
	V100 = lookup_table[j1,k0,l0,:]
	V200 = lookup_table[j2,k0,l0,:]
	V010 = lookup_table[j0,k1,l0,:]
	V110 = lookup_table[j1,k1,l0,:]
	V210 = lookup_table[j2,k1,l0,:]
	V020 = lookup_table[j0,k2,l0,:]
	V120 = lookup_table[j1,k2,l0,:]
	V220 = lookup_table[j2,k2,l0,:]
	
	V001 = lookup_table[j0,k0,l1,:]
	V101 = lookup_table[j1,k0,l1,:]
	V201 = lookup_table[j2,k0,l1,:]
	V011 = lookup_table[j0,k1,l1,:]
	V111 = lookup_table[j1,k1,l1,:]
	V211 = lookup_table[j2,k1,l1,:]
	V021 = lookup_table[j0,k2,l1,:]
	V121 = lookup_table[j1,k2,l1,:]
	V221 = lookup_table[j2,k2,l1,:]
	
	V002 = lookup_table[j0,k0,l2,:]
	V102 = lookup_table[j1,k0,l2,:]
	V202 = lookup_table[j2,k0,l2,:]
	V012 = lookup_table[j0,k1,l2,:]
	V112 = lookup_table[j1,k1,l2,:]
	V212 = lookup_table[j2,k1,l2,:]
	V022 = lookup_table[j0,k2,l2,:]
	V122 = lookup_table[j1,k2,l2,:]
	V222 = lookup_table[j2,k2,l2,:]
	
	
	for i in range(3):
		
		if jp == 'u':
			A = np.array([[1.,j0,j0**2.],[1.,j1,j1**2.],[1.,j2,j2**2.]])
			b = np.array([V000[i],V100[i],V200[i]])
			a = np.linalg.solve(A,b)
			c00 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V010[i],V110[i],V210[i]])
			a = np.linalg.solve(A,b)
			c10 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V020[i],V120[i],V220[i]])
			a = np.linalg.solve(A,b)
			c20 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V001[i],V101[i],V201[i]])
			a = np.linalg.solve(A,b)
			c01 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V011[i],V111[i],V211[i]])
			a = np.linalg.solve(A,b)
			c11 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V021[i],V121[i],V221[i]])
			a = np.linalg.solve(A,b)
			c21 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V002[i],V102[i],V202[i]])
			a = np.linalg.solve(A,b)
			c02 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V012[i],V112[i],V212[i]])
			a = np.linalg.solve(A,b)
			c12 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V022[i],V122[i],V222[i]])
			a = np.linalg.solve(A,b)
			c22 = a[0] + a[1] * j + a[2] * j**2.
		else:
			A = np.array([[1.,j2,j2**2.],[1.,j0,j0**2.],[1.,j1,j1**2.]])
			b = np.array([V200[i],V000[i],V100[i]])
			a = np.linalg.solve(A,b)
			c00 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V210[i],V010[i],V110[i]])
			a = np.linalg.solve(A,b)
			c10 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V220[i],V020[i],V120[i]])
			a = np.linalg.solve(A,b)
			c20 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V201[i],V001[i],V101[i]])
			a = np.linalg.solve(A,b)
			c01 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V211[i],V011[i],V111[i]])
			a = np.linalg.solve(A,b)
			c11 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V221[i],V021[i],V121[i]])
			a = np.linalg.solve(A,b)
			c21 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V202[i],V002[i],V102[i]])
			a = np.linalg.solve(A,b)
			c02 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V212[i],V012[i],V112[i]])
			a = np.linalg.solve(A,b)
			c12 = a[0] + a[1] * j + a[2] * j**2.
			
			b = np.array([V222[i],V022[i],V122[i]])
			a = np.linalg.solve(A,b)
			c22 = a[0] + a[1] * j + a[2] * j**2.
		
		
		
		
		if kp == 'u':
			A = np.array([[1.,k0,k0**2.],[1.,k1,k1**2.],[1.,k2,k2**2.]])
			b = np.array([c00,c10,c20])
			a = np.linalg.solve(A,b)
			c0 = a[0] + a[1] * k + a[2] * k**2.
			
			b = np.array([c01,c11,c21])
			a = np.linalg.solve(A,b)
			c1 = a[0] + a[1] * k + a[2] * k**2.
			
			b = np.array([c02,c12,c22])
			a = np.linalg.solve(A,b)
			c2 = a[0] + a[1] * k + a[2] * k**2.
		else:
			A = np.array([[1.,k2,k2**2.],[1.,k0,k0**2.],[1.,k1,k1**2.]])
			b = np.array([c20,c00,c10])
			a = np.linalg.solve(A,b)
			c0 = a[0] + a[1] * k + a[2] * k**2.
			
			b = np.array([c21,c01,c11])
			a = np.linalg.solve(A,b)
			c1 = a[0] + a[1] * k + a[2] * k**2.
			
			b = np.array([c22,c02,c12])
			a = np.linalg.solve(A,b)
			c2 = a[0] + a[1] * k + a[2] * k**2.
		
		
		if lp == 'u':
			A = np.array([[1.,l0,l0**2.],[1.,l1,l1**2.],[1.,l2,l2**2.]])
			b = np.array([c0,c1,c2])
			a = np.linalg.solve(A,b)
			v[i] = a[0] + a[1] * l + a[2] * l**2.
		else:
			A = np.array([[1.,l2,l2**2.],[1.,l0,l0**2.],[1.,l1,l1**2.]])
			b = np.array([c2,c0,c1])
			a = np.linalg.solve(A,b)
			v[i] = a[0] + a[1] * l + a[2] * l**2.
		
	
	v = v / a1
	
	return v
