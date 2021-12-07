from numpy import pi,sqrt,cos,sin
#a,b,n,c1,c2,phi,w,xp,yp,zp=const
# 1 2 3  4  5   6 7  8  9  10
def nu(a1,a2,b1,b2,n,c1,c2,phi,w,xp,yp,zp,m,p,i,j):
	const1 = (a1,b1,n,c1,c2,phi,w,xp,yp,zp)
	v1 = integrate(const1,m,p)
	const2 = (a2,b2,n,c1,c2,phi,w,xp,yp,zp)
	v2 = integrate(const2,m,p)
	if i==j:
		return v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2]
	else:
		r1 = (xp-a1*cos(phi), yp-a1*sin(phi), zp)
		r2 = (xp-a2*cos(phi), yp-a2*sin(phi), zp)
		t = straight_segment(r1,r2)
		return v2[0]-v1[0]+t[0], v2[1]-v1[1]+t[1], v2[2]-v1[2]+t[2]

def straight_segment(r1,r2):
	r10,r11,r12=r1
	r20,r21,r22=r2
	r1mag = sqrt(r10*r10+r11*r11+r12*r12)
	r2mag = sqrt(r20*r20+r21*r21+r22*r22)
	cross0 = r11*r22-r12*r21
	cross1 = r12*r20-r10*r22
	cross2 = r10*r21-r11*r20
	dot = r10*r20+r11*r21+r12*r22
	k = (r1mag + r2mag) / r1mag / r2mag / (r1mag * r2mag + dot) / 4. / pi
	return k*cross0, k*cross1, k*cross2

def integrate(const,m,p):
	n=const[2]
	w=const[6]
	steps = int(n*m)
	vel = (0.,0.,0.)
	for i in range(steps):
		zeta0 = (float(i) / float(steps)) ** p
		zeta  = (float(i+1) / float(steps)) ** p
		dzeta = zeta - zeta0
		vel = rk4(const,zeta0,dzeta,vel)
	return -w/4./pi*vel[0], -w/4./pi*vel[1], -w/4./pi*vel[2]

def rk4(const,zeta,dzeta,y):
	y0,y1,y2=y
	k10,k11,k12=f(const,zeta)
	k20,k21,k22=f(const,zeta+dzeta/2.)
	k40,k41,k42=f(const,zeta+dzeta)
	h=dzeta/6.
	return y0+h*(k10+4.*k20+k40),y1+h*(k11+4.*k21+k41),y2+h*(k12+4.*k22+k42)

def f(const,zeta):
	den = denom(const,zeta)
	return fx(const,zeta,den), fy(const,zeta,den), fz(const,zeta,den)

def denom(const,zeta):
	c=const[:7]
	xp,yp,zp = const[7:10]
	return ((xp-Xh(c,zeta))**2.+(yp-Yh(c,zeta))**2.+(zp-Zh(c,zeta))**2.)**1.5

def fx(const,zeta,den):
	c=const[:7]
	yp,zp = const[8:10]
	return (Yh_(c,zeta)*(zp-Zh(c,zeta))-Zh_(c)*(yp-Yh(c,zeta)))/den

def fy(const,zeta,den):
	c=const[:7]
	xp,yp,zp=const[7:10]
	return (Zh_(c)*(xp-Xh(c,zeta))-Xh_(c,zeta)*(zp-Zh(c,zeta)))/den

def fz(const,zeta,den):
	c=const[:7]
	xp,yp=const[7:9]
	return (Xh_(c,zeta)*(yp-Yh(c,zeta))-Yh_(c,zeta)*(xp-Xh(c,zeta)))/den

def Xh(c,zeta):
	a,b,n,c1,c2,phi=c[:6]
	z=b*n*zeta
	return a*(1.+(c2-1.)*z/sqrt(z**2.+c1**2.))*cos(2.*pi*n*zeta+phi)

def Yh(c,zeta):
	a,b,n,c1,c2,phi,w=c
	z=b*n*zeta
	return -a*(1.+(c2-1.)*z/sqrt(z**2.+c1**2.))*sin(2.*pi*n*zeta+phi)*w

def Zh(c,zeta):
	b,n=c[1:3]
	return b*n*zeta

def Xh_(c,zeta):
	a,b,n,c1,c2,phi=c[:6]
	z=b*n*zeta
	theta = 2.*pi*n
	d = z**2.+c1**2.
	num = c2 - 1.
	e = sqrt(d)
	return cos(theta*zeta+phi)*a*b*n*num/e*(1.-z**2./d)-a*theta*(1.+z*num/e)*sin(theta*zeta+phi)

def Yh_(c,zeta):
	a,b,n,c1,c2,phi,w=c
	z=b*n*zeta
	theta = 2.*pi*n
	d = z**2.+c1**2.
	num = c2 - 1.
	e = sqrt(d)
	return -w*sin(theta*zeta+phi)*a*b*n*num/e*(1.-z**2./d)-w*a*theta*(1.+z*num/e)*cos(theta*zeta+phi)

def Zh_(c):
	b,n=c[1:3]
	return b*n

def vec_mag(r):
	return sqrt(r[0]**2.+r[1]**2.+r[2]**2.)

