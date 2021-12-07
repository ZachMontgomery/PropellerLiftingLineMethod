from numpy import sin, cos, pi
import numpy as np

def fx(a, b, x, y, z, rl, n, zeta, phi):
	theta = 2*pi*n*zeta + phi
	num = -b*y-rl*a*b*sin(theta)-rl*2*pi*a*z*cos(theta) + rl*2*pi*a*b*n*zeta*cos(theta)
	den = (x-a*cos(theta)) ** 2 + (y + rl*a*sin(theta)) ** 2 + (z - b*n*zeta) ** 2
	vx = num / den ** 1.5
	return vx

def fy(a, b, x, y, z, rl, n, zeta, phi):
	theta = 2*pi*n*zeta + phi
	num = b*x - a*b*cos(theta) + 2*pi*a*z*sin(theta) - 2*pi*a*b*n*zeta*sin(theta)
	den = (x-a*cos(theta)) ** 2 + (y + rl*a*sin(theta)) ** 2 + (z - b*n*zeta) ** 2
	vy = num / den ** 1.5
	return vy

def fz(a, b, x, y, z, rl, n, zeta, phi):
	theta = 2*pi*n*zeta + phi
	num = -rl*a + rl*x*cos(theta) - y*sin(theta)
	den = (x-a*cos(theta)) ** 2 + (y + rl*a*sin(theta)) ** 2 + (z - b*n*zeta) ** 2
	vz = num / den ** 1.5
	return vz

def rk4_v(a, b, x, y, z, n, m, p, rl, phi):
	k = n*m
	nux = 0.
	nuy = 0.
	nuz = 0.
	for i in range(int(k)):
		zeta1 = (i / (k)) ** p
		zeta2 = ((i+1) / (k)) ** p
		dzeta = zeta2 - zeta1
		# calculate x component
		k1 = fx(a, b, x, y, z, rl, n, zeta1, phi)
		k2 = fx(a, b, x, y, z, rl, n, dzeta/2+zeta1, phi)
		#k3 = k2
		k4 = fx(a, b, x, y, z, rl, n, zeta2, phi)
		nux = nux + dzeta/6*(k1+4*k2+k4)
		# calculate y component
		k1 = fy(a, b, x, y, z, rl, n, zeta1, phi)
		k2 = fy(a, b, x, y, z, rl, n, dzeta/2+zeta1, phi)
		#k3 = k2
		k4 = fy(a, b, x, y, z, rl, n, zeta2, phi)
		nuy = nuy + dzeta/6*(k1+4*k2+k4)
		# calculate z component
		k1 = fz(a, b, x, y, z, rl, n, zeta1, phi)
		k2 = fz(a, b, x, y, z, rl, n, dzeta/2+zeta1, phi)
		#k3 = k2
		k4 = fz(a, b, x, y, z, rl, n, zeta2, phi)
		nuz = nuz + dzeta/6*(k1+4*k2+k4)
	v = np.zeros(3)
	v[0] = -rl*nux*n/4./pi
	v[1] = -rl*nuy*n/4./pi
	v[2] = -rl*nuz*a*n/2.
	#v = -rl*np.array([nux*n/4/pi, nuy*n/4/pi, nuz*a*n/2])
	return v

def straight_segment_v(r1, r2):
	r1mag = vec_mag(r1)
	r2mag = vec_mag(r2)
	v = (r1mag + r2mag) * np.cross(r1,r2) / r1mag / r2mag / (r1mag * r2mag + np.dot(r1,r2)) / 4. / pi
	return v

def vec_mag(r):
	mag = (r[0]**2.+r[1]**2.+r[2]**2.)**.5
	return mag

def nu(a1, a2, b, x, y, z, n, m, p, rl, r1, r2, i, j, phi):
	temp = np.zeros(3)
	if i != j:
		temp = -rl*straight_segment_v(r1,r2)
	
	v1 = rk4_v(a1, b, x, y, z, n, m, p, rl, phi)
	v2 = rk4_v(a2, b, x, y, z, n ,m ,p ,rl, phi)
	x = v2 - v1 + temp
	return x


