import numpy as np
import re
import helix as hx
import helix_func as hf
#import lookup_table_reader as ltr

numeric_const_pattern = r"""
     [-+]? # optional sign
     (?:
         (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
         |
         (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
     )
     # followed by optional exponent part if desired
     (?: [Ee] [+-]? \d+ ) ?
     """
rx = re.compile(numeric_const_pattern, re.VERBOSE)
del numeric_const_pattern

def myfunc(f):
	dummy = rx.findall(f.readline())
	return dummy[0]

class propeller:
	
	def __init__(self, filename):
		self.f = open(filename)
		# number of blades
		self.n_blades = int(myfunc(self.f))
		# number of control points per blade
		self.nb = int(myfunc(self.f))
		# total number of control points
		self.cp_total = self.n_blades * self.nb
		# propeller radius
		self.rp = float(myfunc(self.f))
		# hub radius
		self.rh = float(myfunc(self.f))
		# right-handed or left-handed propeller
		self.rl = float(myfunc(self.f))
		# freestream velocity
		self.vinf = float(myfunc(self.f))
		# angular velocity of propeller
		self.w = float(myfunc(self.f)) / 60. * 2. * np.pi
		# number of loops in helices
		self.n = float(myfunc(self.f))
		# avg number of steps per loop
		self.m = float(myfunc(self.f))
		# power clustering factor
		self.p = float(myfunc(self.f))
		# relaxation factor for newtons method
		self.Omega = float(myfunc(self.f))
		# relaxation factor for pitch iteration
		self.Omega1 = float(myfunc(self.f))
		# fluid density
		self.rho = float(myfunc(self.f))
		self.f.close()
		# helical pitch
		self.b = self.vinf / self.w * 2. * np.pi
		if self.b == 0.: self.b = 0.01
		# freestream velocity vector
		self.Vinf = np.array([0., 0., self.vinf])
		# angular velocity vector
		if self.rl == -1.:
			self.W = np.array([0., 0., -self.w])
			print('Right-handed prop')
		else:
			self.W = np.array([0., 0., self.w])
			print('Left-handed prop')
		# advance ratio
		self.J = np.pi * self.vinf / self.w / self.rp
		
		# initialize node and control point matricies to 0, and chord and twist matrices
		self.node = np.zeros((self.cp_total,2,3))
		self.c = np.zeros((self.cp_total,2))
		self.cp = np.zeros((self.cp_total,3))
		self.twist = np.zeros((self.cp_total))
		# initialize node1 and chord on blade one
		self.node[0,0,0] = self.rh
		self.node[0,0,1] = 0.
		self.c[0,0] = .1
		# mirror node1 and chord on blade one to the other blades
		angle_spacing = 2. * np.pi / float(self.n_blades)
		for i in range(1,self.n_blades):
			phi = float(i) * angle_spacing
			self.node[i*self.nb,0,0] = self.rh * np.cos(phi)
			self.node[i*self.nb,0,1] = self.rh * np.sin(phi)
			self.c[i*self.nb,0] = self.c[0,0]
		# loop through for the rest of the points on blade one
		for i in range(self.nb):
			# calculate control point location on first blade using cosine clustering and calculate twist
			self.cp[i,0] = (self.rp - self.rh) / 2. * (1. - np.cos(float(i+1)*np.pi/self.nb - np.pi / 2. / self.nb)) + self.rh
			self.cp[i,1] = 0.
			self.twist[i] = .5 * np.pi + 10. * np.pi / 180.
			# calculate node2 point location on first blade using cosine clustering
			self.node[i,1,0] = (self.rp - self.rh) / 2. * (1. - np.cos(float(i+1) * np.pi / self.nb)) + self.rh
			self.node[i,1,1] = 0.
			self.c[i,1] = .1
			# calculate node1 point location on first blade using cosine clustering
			if i != self.nb - 1:
				self.node[i+1,0,:] = self.node[i,1,:]
				self.c[i+1,0] = self.c[i,1]
			# mirror these control and node points onto the other blades
			for j in range(1,self.n_blades):
				phi = float(j) * angle_spacing
				# control point
				r = self.cp[i,0]
				self.cp[i+j*self.nb,0] = r * np.cos(phi)
				self.cp[i+j*self.nb,1] = r * np.sin(phi)
				self.twist[i+j*self.nb] = self.twist[i]
				# node1
				r = self.node[i,0,0]
				self.node[i+j*self.nb,0,0] = r * np.cos(phi)
				self.node[i+j*self.nb,0,1] = r * np.sin(phi)
				self.c[i+j*self.nb,0] = self.c[i,0]
				# node 2
				r = self.node[i,1,0]
				self.node[i+j*self.nb,1,0] = r * np.cos(phi)
				self.node[i+j*self.nb,1,1] = r * np.sin(phi)
				self.c[i+j*self.nb,1] = self.c[i,1]
		#del r, phi, angle_spacing
		
		# calculate ri and Vpi vectors
		self.Vp = np.zeros((self.cp_total,3))
		self.r = np.copy(self.cp[:,:])
		for i in range(self.cp_total):
			self.Vp[i,:] = self.Vinf + np.cross(self.r[i,:],self.W)
		
		# calculate nu_ji
		#self.renew_nu()
		
		# initialize the unit vectors
		self.un = np.zeros((self.cp_total,3))
		self.ua = np.zeros((self.cp_total,3))
		negk = np.array([0.,0.,-1.])
		k = -negk[:]
		# loop through the control points to calculate the unit vectors
		for i in range(self.cp_total):
			temp = np.cross(self.r[i],self.W)
			temp = temp[:] / hx.vec_mag(temp)
			self.un[i,:] = np.cos(self.twist[i]) * negk[:] + np.sin(self.twist[i]) * temp[:]
			self.ua[i,:] = np.cos(self.twist[i]) * temp[:] + np.sin(self.twist[i]) * k[:]
		
		# calculate dA
		self.dA = np.zeros((self.cp_total))
		for i in range(self.cp_total):
			self.dA[i] = (self.c[i,0] + self.c[i,1]) / 2. * (hx.vec_mag(self.node[i,1,:]) - hx.vec_mag(self.node[i,0,:]))
		
		# calculate dl
		self.dl = np.zeros((self.cp_total,3))
		for i in range(self.cp_total):
			self.dl[i,:] = self.rl*(self.node[i,0,:] - self.node[i,1,:])
	
	def first_guess(self):
		# right side vector
		b = np.zeros((self.cp_total))
		for i in range(self.cp_total):
			b[i] = hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.Vp[i,:],self.un[i,:]) * self.dA[i]
		
		# calculate A matrix
		A = np.zeros((self.cp_total, self.cp_total))
		for i in range(self.cp_total):
			for j in range(self.cp_total):
				if i == j:
					A[i,j] = 2. * hx.vec_mag(np.cross(self.Vp[i,:],self.dl[i,:])) - hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.nu[i,j,:],self.un[i,:]) * self.dA[i]
				else:
					A[i,j] = - hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.nu[i,j,:],self.un[i,:]) * self.dA[i]
		
		# invert A matrix
		#Ainv = np.linalg.inv(A)
		
		# create first guess for the gamma vector
		self.gamma = np.linalg.solve(A,b)#Ainv.dot(b)
	
	def calc_V(self):
		self.V = np.zeros((self.cp_total,3))
		self.aoa = np.zeros((self.cp_total))
		for i in range(self.cp_total):
			temp = np.zeros((3))
			for j in range(self.cp_total):
				temp += self.gamma[j] * self.nu[i,j,:]
			self.V[i,:] = self.Vp[i,:] + temp
			self.aoa[i] = np.arctan(np.dot(self.V[i,:],self.un[i,:])/np.dot(self.V[i,:],self.ua[i,:]))
			if self.aoa[i] < 0.:
				print('Angle of attack is {} at control point {}'.format(self.aoa[i], i))
	
	def calc_Cl(self):
		self.Cl = np.zeros((self.cp_total))
		self.Cla = np.zeros((self.cp_total))
		for i in range(self.cp_total):
			if self.aoa[i] <= .25:
				self.Cl[i] = 2. * np.pi * self.aoa[i]
				self.Cla[i] = 2. * np.pi
			else:
				print('Section {} stalled with {} degrees aoa.'.format(i,self.aoa[i]*180./np.pi))
				self.Cl[i] = np.pi / 2. * np.cos(self.aoa[i]) / np.cos(.25)
				self.Cla[i] = -np.pi / 2. * np.sin(self.aoa[i]) / np.cos(.25)
	
	def calc_R(self):
		self.R = np.zeros((self.cp_total))
		for i in range(self.cp_total):
			self.R[i] = 2.*self.gamma[i]*hx.vec_mag(np.cross(self.V[i,:],self.dl[i,:]))-hx.vec_mag(self.V[i,:])**2.*self.Cl[i]*self.dA[i]
	
	def calc_Jacobian(self):
		# initialize vectors
		w = np.zeros((self.cp_total,3))
		Vn = np.zeros((self.cp_total))
		Va = np.zeros((self.cp_total))
		self.Jacobian = np.zeros((self.cp_total,self.cp_total))
		for i in range(self.cp_total):
			w[i,:] = np.cross(self.V[i,:],self.dl[i,:])
			Vn[i] = np.dot(self.V[i,:],self.un[i,:])
			Va[i] = np.dot(self.V[i,:],self.ua[i,:])
			for j in range(self.cp_total):
				term1 = 2. * self.gamma[i] * np.dot(w[i,:],np.cross(self.nu[i,j,:],self.dl[i,:])) / hx.vec_mag(w[i,:])
				
				term2 = 2. * self.Cl[i] * np.dot(self.V[i,:],self.nu[i,j,:]) * self.dA[i]
				
				term3 = hx.vec_mag(self.V[i,:]) ** 2. * self.Cla[i] * (Va[i] * np.dot(self.nu[i,j,:],self.un[i,:]) - Vn[i] * np.dot(self.nu[i,j,:],self.ua[i,:])) / (Vn[i] ** 2. + Va[i] ** 2.) * self.dA[i]
				
				term4 = 2. * hx.vec_mag(w[i,:])
				if i != j:
					self.Jacobian[i,j] = term1 - term2 - term3
				else:
					self.Jacobian[i,j] = term1 - term2 - term3 + term4
	
	def update_gamma(self, rf):
		dgamma = np.zeros((self.cp_total))
		dgamma = np.linalg.solve(self.Jacobian,-self.R)
		self.gamma = self.gamma + rf * dgamma
	
	def calc_gamma(self, tol, max_it):
		self.calc_V()
		self.calc_Cl()
		self.calc_R()
		i = 0
		b = 0
		e = max(abs(self.R))
		print(e, i)
		temp = self.Omega
		while e > tol:
			self.calc_Jacobian()
			self.update_gamma(temp)
			self.calc_V()
			self.calc_Cl()
			self.calc_R()
			e = max(abs(self.R))
			i += 1
			print(e, i)
			print()
			if i >= max_it:
				print('Newtons method failed to converge in {} iterations!'.format(i))
				#input()
				temp -= .1
				i = 0
				b += 1
				if b == 5:
					break
				#break
		print('final gammas ',self.gamma)
		#input()
	
	def calc_forces_and_moments(self):
		self.dF = np.zeros((self.cp_total, 3))
		self.total_force = np.zeros(3)
		self.dM = np.zeros((self.cp_total, 3))
		self.total_moment = np.zeros(3)
		self.us = np.zeros((self.cp_total, 3))
		self.Cm = np.zeros(self.cp_total)
		for i in range(self.cp_total):
			self.dF[i,:] = self.rho * self.gamma[i] * np.cross(self.V[i,:],self.dl[i,:])
			self.total_force += self.dF[i,:]
			
			self.us[i,:] = np.cross(self.ua[i,:],self.un[i,:])
			self.Cm[i] = -.1
			self.dM[i,:] = -1./8.*self.rho*hx.vec_mag(self.V[i,:])**2.*self.Cm[i]*(self.c[i,0]+self.c[i,1])**2.*(hx.vec_mag(self.node[i,1,:])-hx.vec_mag(self.node[i,0,:]))*self.us[i,:]
			self.total_moment += np.cross(self.r[i,:],self.dF[i,:]) + self.dM[i,:]
		print('\n\nThe total force is  {}'.format(self.total_force))
		print('The total moment is {}\n\n'.format(self.total_moment))
	
	def calc_coef(self):
		self.CT = -self.total_force[2] / self.rho / (2.*self.rp)**4. / (self.w/2./np.pi)**2.
		self.CL = -self.rl * self.total_moment[2] / abs(self.total_moment[2]) * hx.vec_mag(self.total_moment) / self.rho / (2.*self.rp)**5. / (self.w/2./np.pi)**2.
		self.CP = self.CL * 2. * np.pi
		print('Coefficient of thrust is {}'.format(self.CT))
		print('Coefficient of torque is {}'.format(self.CL))
		print('Coefficient of power is  {}\n\n'.format(self.CP))
	
	def pll(self, tol, max_it):
		self.first_guess()
		self.calc_gamma(tol, max_it)
		self.calc_forces_and_moments()
		self.calc_coef()
	
	def renew_nu(self):
		# calculate nu_ji
		self.nu = np.zeros((self.cp_total,self.cp_total, 3))
		for j in range(self.cp_total):
			a1 = (self.node[j,0,0] ** 2. + self.node[j,0,1] ** 2.) ** .5
			a2 = (self.node[j,1,0] ** 2. + self.node[j,1,1] ** 2.) ** .5
			phi = np.arctan2(self.cp[j,1],self.cp[j,0])
			#print('phi is {}'.format(phi))
			for i in range(self.cp_total):
				#v1 = np.zeros(3)
				#v2 = np.zeros(3)
				#vh1 = np.zeros(3)
				#vh2 = np.zeros(3)
				#st = np.zeros(3)
				#if phi != 0.:
					#r = (self.cp[i,0]**2. + self.cp[i,1]**2.)**.5
					#theta = np.arctan2(self.cp[i,1],self.cp[i,0]) - phi
					#cpx = r * np.cos(theta)
					#cpy = r * np.sin(theta)
					#v1 = ltr.read_table(a1, self.b, cpx, cpy)
					#vh1[0] = v1[0] * np.cos(phi) + v1[1] * np.sin(phi)
					#vh1[1] = v1[1] * np.cos(phi) - v1[0] * np.sin(phi)
					#vh1[2] = v1[2]
					#v2 = ltr.read_table(a2, self.b, cpx, cpy)
					#vh2[0] = v2[0] * np.cos(phi) + v2[1] * np.sin(phi)
					#vh2[1] = v2[1] * np.cos(phi) - v2[0] * np.sin(phi)
					#vh2[2] = v2[2]
				#else:
					#cpx = self.cp[i,0]
					#cpy = self.cp[i,1]
					#vh1 = ltr.read_table_quadratic(a1, self.b, cpx, cpy)
					#vh2 = ltr.read_table_quadratic(a2, self.b, cpx, cpy)
				
				#if i != j:
					#r1 = self.cp[i,:] - self.node[j,0,:]
					#r2 = self.cp[i,:] - self.node[j,1,:]
					#st = hf.straight_segment(st,r1,r2)
				
				
				#v = vh2 - vh1 - self.rl * st
				#print(v,'lookup table')
				
				w = hf.calc_nu(self.nu[i,j,:],a1,a2,self.b,self.rl,self.cp[i,:],self.n,self.m,self.p,phi,'r',i==j)
				#print(w, 'python helix func')
				#error = abs((w - v) / w) * 100.
				#print(error,'% error')
				
				#input()
				#print('Finished horseshoe vortex {0} on control point {1}'.format(j+1,i+1))
			print('Finished calculating velocity induced by horseshoe vortex {}'.format(j+1))
	
	def renew_Vp(self):
		self.Vp = np.zeros((self.cp_total,3))
		for i in range(self.cp_total):
			self.Vp[i,:] = self.Vinf + np.cross(self.r[i,:],self.W)
	
