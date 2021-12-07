import numpy as np

def chord(rp, r):
	c = .075 * 2. * rp * (1. - (r / rp) ** 2.) ** .5
	return c

def twist(r, global_pitch):
	t = np.arctan((1.8 - 2.*np.pi*r*np.tan(-2.1*np.pi/180.))/(2.*np.pi*r + 1.8*np.tan(-2.1*np.pi/180.))) + global_pitch
	return t

def Cl(aoa):
	
	if aoa <= .25:
		cl = 2. * np.pi * aoa
		cla = 2. * np.pi
	else:
		cl = np.pi / 2. * np.cos(aoa) / np.cos(.25)
		cla = -np.pi / 2. * np.sin(aoa) / np.cos(.25)  #2. * np.pi 
	
	#~ if aoa <= .25:
		#~ cl = 2. * np.pi * aoa
		#~ cla = 2. * np.pi
	#~ else:
		#~ cl = 0.
		#~ cla = 0.
	
	#~ cl = 2. * np.pi * aoa
	#~ cla = 2. * np.pi
	
	#~ k1 = 2. / np.pi
	#~ k2 = -0.36549726951209593
	#~ k3 = 5.0753455211144347
	#~ if aoa <= .25:
		#~ cl = 2. * np.pi * aoa
		#~ cla = 2. * np.pi
	#~ elif aoa <= 1.0005006021307872:
		#~ cl  = k1 * np.pi * np.sin( aoa * k3 + k2 )
		#~ cla = k3 * k1 * np.pi * np.cos( aoa * k3 + k2 )
	#~ else:
		#~ cl  = k1 * np.pi * np.sin( 1.0005006021307872 * k3 + k2 )
		#~ cla = 0.
	
	#~ k = 5.263157894736842
	#~ cl = .5*np.pi*np.sin((aoa+2.1*np.pi/180.)*k)
	#~ cla = .5*np.pi*k*np.cos((aoa+2.1*np.pi/180.)*k)
	
	#~ cl = -aoa * 2.5 * (aoa - np.pi/2.)
	#~ cla = -2. * 2.5 * aoa + np.pi/2. * 2.5
	
	return cl, cla


def Cd(aoa):
	if aoa <= .25:
		cd = 0.244 * aoa **2. + 0.006
	elif aoa <= 0.3:
		cd = 16.6944 * aoa ** 2. - 1.0234
	else:
		cd = np.pi / 2. * np.sin(aoa) / np.cos(.25)
	return cd


