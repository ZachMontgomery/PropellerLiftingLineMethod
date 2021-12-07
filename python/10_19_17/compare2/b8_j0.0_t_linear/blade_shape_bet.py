from numpy import pi
from numpy import sin, cos, tan, arctan, arccos, exp


# chord function
def chord(dp, zeta):
	cb = .075*dp*(1.-zeta**2.)**.5
	return cb

# function
def func(cb_, zeta, CL, k, Bt, ei, einf):
	f = cb_/8./zeta*CL - arccos(exp(-k*(1.-zeta)/2./sin(Bt)))*tan(ei)*sin(einf+ei)
	# print(' f is {}'.format(f))
	return f

# derivative of function
def dfunc(cb_, zeta, CL_ei, k, Bt, ei, einf):
	part1 = -arccos(exp(-k*(1.-zeta)/2./sin(Bt)))*tan(ei)*cos(einf+ei)
	part2 = -arccos(exp(-k*(1.-zeta)/2./sin(Bt)))/cos(ei)**2.*sin(einf+ei)
	part3 = cb_/8./zeta*CL_ei
	df = part1 + part2 + part3
	# print('			part1 is {}\n			part2 is {}'.format(part1,part2))
	if df == 0.:
		print('Manually set df to 1')
		df = 1.
	return df

# coefficient of lift
def Cl(ab):
	# if ab > .25:
		# lift = pi/2.*cos(ab)/cos(.25)
	# else:
		# lift = 2.*pi*ab
	
	lift = -ab * 2.5 * (ab - pi / 2.)
	
	return lift

# derivative of lift coefficient
def Cl_ei(ab):
	# if ab > .25:
		# l = pi/2/cos(.25)*sin(ab)
	# else:
		# l = -2.*pi
	
	l = 5. * ab - 5. / 4. * pi
	
	return l

# calculate angle of attack
def aoa(B, einf, ei):
	ab = B - einf - ei
	return ab

# newtons method to find ei
def newtons_method(xo, tol, cb_, zeta, k, Bt, einf, B):
	ea = 100.
	j = 0.
	while ea > tol:
		j += 1.
		ab = aoa(B, einf, xo)
		CL = Cl(ab)
		CL_ei = Cl_ei(ab)
		
		f = func(cb_, zeta, CL, k, Bt, xo, einf)
		
		df = dfunc(cb_, zeta, CL_ei, k, Bt, xo, einf)
		xn = xo - (f / df)
		ea = abs((xn - xo) / xo)
		xo = xn
		if j > 200.:
			print('\n\n\nToo many iterations in Newtons method!\n\n\n')
			break
	# input()
	return xn

# drag coefficient
def Cd(ab):
	if ab <= .25:
		d = .224*ab**2.+.006
	elif ab <= .3:
		d = 16.6944*ab**2.-1.0234
	else:
		d = pi/2.*sin(ab)/cos(.25)
	# d = 0.
	return d

# ctp
def Ctp(zeta, einf, ei, ap, cb_, cb, dl):
	CL = Cl(ap)
	CD = Cd(ap)
	Ct = pi**2./4.*zeta**2.*cb_*(cos(ei)/cos(einf))**2.*(CL*cos(einf+ei)-CD*sin(einf+ei))
	
	rho = .0023769
	w = 2400. / 60. * 2. * pi
	k = 30.
	
	z = -rho*w**2./2.*zeta**2.*cb*(cos(ei)/cos(einf))**2.*(CL*cos(einf+ei) - CD*sin(einf+ei)) * dl
	y = rho*w**2./2.*zeta**2.*cb*(cos(ei)/cos(einf))**2.*(CD*cos(einf+ei) + CL*sin(einf+ei)) * dl
	
	return Ct, z, y

# Clp
def Clp(zeta, einf, ei, ap, cb_):
	CL = Cl(ap)
	CD = Cd(ap)
	Cll = pi**2./8.*zeta**3.*cb_*(cos(ei)/cos(einf))**2.*(CD*cos(einf+ei)+CL*sin(einf+ei))
	return Cll

# trap
def trap(x, f):
	j = x.size
	area = 0.
	for i in range(j-1):
		area += (x[i+1]-x[i]) * (f[i+1]+f[i])/2.
	return area

def twist(zeta, kc):
	# t = arctan((kc - pi*zeta*tan(-2.1*pi/180.))/(pi*zeta + kc*tan(-2.1*pi/180.)))
	t = ((1. - zeta) * kc + 10.) * pi / 180.
	# t = 10. * pi / 180.
	return t
