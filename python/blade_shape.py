import numpy as np

def chord(rp, r):
	c = .075 * 2. * rp * (1. - (r / rp) ** 2.) ** .5
	return c

def twist(r):
	t = np.arctan((1.8 - 2.*np.pi*r*np.tan(-2.1*np.pi/180.))/(2.*np.pi*r + 1.8*np.tan(-2.1*np.pi/180.)))
	return t

def Cl(aoa):
	if aoa <= .25:
		cl = 2. * np.pi * aoa
		cla = 2. * np.pi
	else:
		print('Stalled with {:5.2f} degrees angle of attack'.format(aoa * 180. / np.pi))
		cl = np.pi / 2. * np.cos(aoa) / np.cos(.25)
		cla = -np.pi / 2. * np.sin(aoa) / np.cos(.25)
	return cl, cla




