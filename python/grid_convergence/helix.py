from numpy import sin, cos, pi
import numpy as np

def fx(zeta,denom,cx1,cx2,cx3,cx4,s,c):
	dx = (cx1 + cx2 * s + (cx3 + cx4 * zeta) * c) / denom
	return dx

def fy(zeta,denom,cy1,cy2,cy3,cy4,s,c):
	dy = (cy1 + cy2 * c + (cy3 + cy4 * zeta) * s) / denom
	return dy

def fz(denom,cz1,cz2,cz3,s,c):
	dz = (cz1 + cz2 * c + cz3 * s) / denom
	return dz

def rk4_v(a, b, x, y, z, n, m, p, rl, phi):
	steps = int(n*m)
	k1 = np.zeros(3)
	k2 = np.zeros(3)
	k4 = np.zeros(3)
	dv_o = np.zeros(3)
	v = np.zeros(3)
	# calculate constants
	c0 = 2. * pi * n
	c1 = n / 4. / pi
	c2 = a * n / 2.
	# x function constants
	cx1 = -b * y
	cx2 = -rl * a * b
	cx3 = -rl * 2. * pi * a * z
	cx4 = rl * 2. * pi * a * b * n
	# y function constants
	cy1 = b * x
	cy2 = -a * b
	cy3 = 2. * pi * a * z
	cy4 = -2. * pi * a * b * n
	# z function constants
	cz1 = -rl * a
	cz2 = rl * x
	cz3 = -y
	# denom constants
	d1 = b * n
	d2 = rl * a
	# calculate zetas
	dzeta = 1. / n / m
	zeta = np.zeros(steps + 1)
	for i in range(steps + 1):
		zeta[i] = (float(i) * dzeta) ** p
	
	theta = phi
	s = sin(theta)
	c = cos(theta)
	denom = ((x-a*c)**2. + (y+d2*s)**2. + z**2.)**1.5
	
	dv_o[0] = fx(0.,denom,cx1,cx2,cx3,cx4,s,c)
	dv_o[1] = fy(0.,denom,cy1,cy2,cy3,cy4,s,c)
	dv_o[2] = fz(denom,cz1,cz2,cz3,s,c)
	
	
	for i in range(1,steps+1):
		dzeta = zeta[i] - zeta[i-1]
		# calculate k1
		k1 = dv_o
		# calculate k2
		zeta_temp = zeta[i-1] + dzeta / 2.
		theta = c0 * zeta_temp + phi
		c = cos(theta)
		s = sin(theta)
		denom = ((x-a*c) ** 2. + (y + d2*s) ** 2. + (z - d1*zeta_temp) ** 2.)**1.5
		k2[0] = fx(zeta_temp,denom,cx1,cx2,cx3,cx4,s,c)
		k2[1] = fy(zeta_temp,denom,cy1,cy2,cy3,cy4,s,c)
		k2[2] = fz(denom,cz1,cz2,cz3,s,c)
		# calculate k4
		theta = c0 * zeta[i] + phi
		s = sin(theta)
		c = cos(theta)
		denom = ((x-a*c) ** 2. + (y + d2*s) ** 2. + (z - d1*zeta[i]) ** 2.)**1.5
		k4[0] = fx(zeta[i],denom,cx1,cx2,cx3,cx4,s,c)
		k4[1] = fy(zeta[i],denom,cy1,cy2,cy3,cy4,s,c)
		k4[2] = fz(denom,cz1,cz2,cz3,s,c)
		# add on this step
		v += dzeta / 6. * (k1 + 4. * k2 + k4)
		# save values for faster run times
		dv_o = k4
	v[0] *= -rl * c1
	v[1] *= -rl * c1
	v[2] *= -rl * c2
	return v

def straight_segment_v(r1, r2):
	r1mag = vec_mag(r1)
	r2mag = vec_mag(r2)
	v = (r1mag + r2mag) * np.cross(r1,r2) / r1mag / r2mag / (r1mag * r2mag + np.dot(r1,r2)) / 4. / pi
	return v

def vec_mag(r):
	mag = np.sqrt(r[0]**2.+r[1]**2.+r[2]**2.)
	return mag

def nu(a1, a2, b, control, n, m, p, rl, ij, phi):
	temp = np.zeros(3)
	x = control[0]
	y = control[1]
	z = control[2]
	
	if ij:
		j1 = np.zeros(3)
		j2 = np.zeros(3)
		r1 = np.zeros(3)
		r2 = np.zeros(3)
		
		j1[0] = a1 * cos(phi)
		j1[1] = a1 * sin(phi)
		
		j2[0] = a2 * cos(phi)
		j2[1] = a2 * sin(phi)
		
		r1 = control - j1
		r2 = control - j2
		
		temp = -rl*straight_segment_v(r1,r2)
	
	v1 = rk4_v(a1, b, x, y, z, n, m, p, rl, phi)
	v2 = rk4_v(a2, b, x, y, z, n ,m ,p ,rl, phi)
	x = v2 - v1 + temp
	return x


