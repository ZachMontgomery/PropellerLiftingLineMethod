# blade theory program
from numpy import pi
from numpy import sin, cos, tan, arctan, arccos, exp
import numpy as np
import blade_shape_bet as blade_shape
import matplotlib.pyplot as plt
import json

def blade_element_theory(J,res):
	f = json.load(open('propeller.json'))
	
	dp = f['propeller']['tip radius'] * 2.
	dh = f['propeller']['hub radius'] * 2.
	rp = f['propeller']['tip radius']
	rh = f['propeller']['hub radius']
	tanal0 = tan(-2.1*pi/180.)
	k = f['propeller']['blades']['number']
	Kc = f['propeller']['pitch'] / dp
	temp = pi*(Kc-pi*tanal0)/(pi+Kc*tanal0)
	Bt = blade_shape.twist(1., Kc)
	if Bt <= 1.e-10:
		print('old Bt is',Bt*180./pi)
		input()
		# Bt = Bt + pi/2
		# print('new Bt is',Bt*180./pi)
		# input()
	
	
	c = np.zeros(res)#, dtype = np.float128)
	r = np.zeros(res)#, dtype = np.float128)
	K = np.zeros(res)#, dtype = np.float128)
	zeta = np.zeros(res)#, dtype = np.float128)
	beta = np.zeros(res)#, dtype = np.float128)
	e_inf = np.zeros(res)#, dtype = np.float128)
	ei = np.zeros(res)#, dtype = np.float128)
	ap = np.zeros(res)#, dtype = np.float128)
	Ct = np.zeros(res)#, dtype = np.float128)
	Cl = np.zeros(res)#, dtype = np.float128)
	Cp = np.zeros(res)#, dtype = np.float128)
	z = np.zeros(res)#, dtype = np.float128)
	y = np.zeros(res)#, dtype = np.float128)
	va = np.zeros(res)
	vt = np.zeros(res)
	dl = np.zeros(res)
	node1 = np.zeros(res)
	node2 = np.zeros(res)
	gamma = np.zeros(res)
	
	rho = f['scenario']['density']
	w = f['propeller']['RPMs'] / 60.
	redim = rho * w ** 2. * dp ** 4.
	
	
	node1[0] = rh
	
	
	for i in range(res):
		#zeta[i] = (i / (res-1) * (rp - rh) + rh) / rp
		zeta[i] = ((rp - rh) / 2. * (1. - cos(pi / float(res) * (float(i+1) - .5))) + rh) / rp
		
		node2[i] = ((rp - rh) / 2. * (1. - cos(pi / float(res) * float(i+1))) + rh) / rp
		if i != res - 1:
			node1[i+1] = node2[i]
		
		dl[i] = node2[i] - node1[i]
		
		if zeta[i] >= 1.:
			Ct[i] = 0.
			Cl[i] = 0.
		else:
			c[i] = blade_shape.chord(dp,zeta[i])
			r[i] = zeta[i]*rp
			K[i] = pi*zeta[i]*(Kc-pi*zeta[i]*tanal0)/(pi*zeta[i]+Kc*tanal0)
			beta[i] = blade_shape.twist(zeta[i], Kc)
			e_inf[i] = arctan(J/pi/zeta[i])
			ei[i] = blade_shape.newtons_method(20.*pi/180., 1E-15, k*c[i]/dp, zeta[i], k, Bt, e_inf[i], beta[i])
			ap[i] = blade_shape.aoa(beta[i], e_inf[i], ei[i])
			Ct[i], z[i], y[i] = blade_shape.Ctp(zeta[i], e_inf[i], ei[i], ap[i], k*c[i]/dp, c[i], dl[i])
			Cl[i] = blade_shape.Clp(zeta[i], e_inf[i], ei[i], ap[i], k*c[i]/dp)
			va[i] = 2. * pi * w * r[i] * sin(ei[i]) * cos(ei[i] + e_inf[i]) / cos(e_inf[i])
			vt[i] = 2. * pi * w * r[i] * sin(ei[i]) * sin(ei[i] + e_inf[i]) / cos(e_inf[i])
			
			gamma[i] = blade_shape.calc_gamma(w*2.*np.pi, zeta[i], c[i], blade_shape.Cl(ap[i]), ei[i], e_inf[i])
			
	
	Cp = Cl * 2. * pi
	
	CT = blade_shape.trap(zeta, Ct)
	
	CP = blade_shape.trap(zeta, Cp)
	
	eta_p = CT*J/CP
	
	CL = blade_shape.trap(zeta, Cl)
	
	Tz = 0.
	Ty = 0.
	for j in range(int(k)):
		for i in range(res):
			Tz += z[i]
			Ty += y[i]
	CT_p = -Tz/redim
	
	#plt.figure()
	#plt.plot(zeta,Ct,'-g')
	#plt.plot(zeta,Ct,'sb')
	#plt.show()
	
	# return zeta, y, z, CT, CL, beta, CT_p
	return gamma

