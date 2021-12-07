import numpy as np
import re
import helix as hx
import helix_func as hf
import blade_shape as bs
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
		# manual pitch
		self.manual_pitch = float(myfunc(self.f)) * np.pi / 180.
		self.f.close()
		# advance ratio
		self.J = np.pi * self.vinf / self.w / self.rp
		# pitch
		self.update_pitch()
		# freestream velocity vector
		self.Vinf = np.array([0., 0., self.vinf])
		# angular velocity vector
		if self.rl == -1.:
			self.W = np.array([0., 0., -self.w])
			print('Right-handed prop')
		else:
			self.W = np.array([0., 0., self.w])
			print('Left-handed prop')
		
		# initialize node and control point matricies to 0, and chord and twist matrices
		self.node = np.zeros((self.cp_total,2,3))
		self.c = np.zeros((self.cp_total,2))
		self.cp = np.zeros((self.cp_total,3))
		self.twist = np.zeros((self.cp_total))
		# initialize node1 and chord on blade one
		self.node[0,0,0] = self.rh
		self.node[0,0,1] = 0.
		self.c[0,0] = bs.chord(self.rp,self.rh)
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
			self.twist[i] = bs.twist(self.cp[i,0], self.manual_pitch)
			# calculate node2 point location on first blade using cosine clustering
			self.node[i,1,0] = (self.rp - self.rh) / 2. * (1. - np.cos(float(i+1) * np.pi / self.nb)) + self.rh
			self.node[i,1,1] = 0.
			self.c[i,1] = bs.chord(self.rp,self.node[i,1,0])
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
		del r, phi, angle_spacing
		
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
	
	def update_pitch(self, **arg):
		
		# Guess for CT
		if 'Ct' in arg:
			self.Ct_guess = arg['Ct']
		else:
			self.Ct_guess = float(input('\n\nAdvance ratio is {}, Enter guess for CT: '.format(self.J)))
		# Vg
		if self.vinf == 0.:
			self.Vg = 2. * self.w * self.rp * np.sqrt(self.Ct_guess / np.pi ** 3.)
		else:
			self.Vg = (np.sqrt(self.vinf ** 2. + 16. * self.Ct_guess * (self.w * self.rp) ** 2. / np.pi ** 3.) - self.vinf) / 2.
		
		# helical pitch
		self.b = np.ones([self.cp_total,2]) * (self.vinf + self.Vg) / self.w * 2. * np.pi
	
	def update_advance_ratio(self, **arg):
		if 'J' in arg:
			self.J = arg['J']
		else:
			self.J = float(input('\n\nEnter value for Advance Ratio: '))
		self.vinf = self.J * self.w * self.rp / np.pi
		self.Vinf[2] = self.vinf
		self.renew_Vp()
		if 'Ct' in arg:
			self.update_pitch(Ct = arg['Ct'])
		else:
			self.update_pitch()
	
	def first_guess(self):
		
		
		# first guess without small angle approximations
		
		b = np.zeros(self.cp_total)
		A = np.zeros((self.cp_total,self.cp_total))
		for i in range(self.cp_total):
			b[i] = self.Vg
			for j in range(self.cp_total):
				A[i,j] = self.nu[i,j,2]
		
		#~ # right side vector
		#~ b = np.zeros((self.cp_total))
		#~ for i in range(self.cp_total):
			#~ b[i] = hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.Vp[i,:],self.un[i,:]) * self.dA[i]
		
		#~ # calculate A matrix
		#~ A = np.zeros((self.cp_total, self.cp_total))
		#~ for i in range(self.cp_total):
			#~ for j in range(self.cp_total):
				#~ if i == j:
					#~ A[i,j] = 2. * hx.vec_mag(np.cross(self.Vp[i,:],self.dl[i,:])) - hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.nu[i,j,:],self.un[i,:]) * self.dA[i]
				#~ else:
					#~ A[i,j] = - hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.nu[i,j,:],self.un[i,:]) * self.dA[i]
		
		self.gamma = np.linalg.solve(A,b)
	
	def calc_V(self):
		self.V = np.zeros((self.cp_total,3))
		self.aoa = np.zeros((self.cp_total))
		for i in range(self.cp_total):
			temp = np.zeros((3))
			for j in range(self.cp_total):
				temp += self.gamma[j] * self.nu[i,j,:]
			self.V[i,:] = self.Vp[i,:] + temp
			self.aoa[i] = np.arctan(np.dot(self.V[i,:],self.un[i,:])/np.dot(self.V[i,:],self.ua[i,:]))
			#if self.aoa[i] < 0.:
				#print('Angle of attack is {} at control point {}'.format(self.aoa[i], i))
	
	def calc_Cl(self):
		self.Cl = np.zeros((self.cp_total))
		self.Cla = np.zeros((self.cp_total))
		for i in range(self.cp_total):
			self.Cl[i], self.Cla[i] = bs.Cl(self.aoa[i])
	
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
		en = e
		print('\n\nBegin Newton\'s solver\n\n')
		print('Iteration = {:<3d}  Max residual = {:6.1e}'.format(i,e))
		temp = self.Omega
		while e > tol:
			eo = en
			self.calc_Jacobian()
			self.update_gamma(temp)
			self.calc_V()
			self.calc_Cl()
			self.calc_R()
			e = max(abs(self.R))
			en = e
			i += 1
			if en > eo:
				temp -= .01
				print('Relaxed Newton\'s Method further to promote convergence. Relaxation factor is {}'.format(temp))
				if temp <= 1.e-10:
					input()
					break
			
			print('Iteration = {:<3d}  Max residual = {:6.1e}'.format(i,e))
			
			if i >= max_it:
				print('Newtons method failed to converge in {} iterations!'.format(i))
				#temp -= .0005
				i = 0
				b += 1
				if b == 5:
					input()
					break
	
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
		self.CT = -self.total_force[2] / abs(self.total_force[2]) * hx.vec_mag(self.total_force)  / self.rho / (2.*self.rp)**4. / (self.w/2./np.pi)**2.
		self.CL = -self.rl * self.total_moment[2] / abs(self.total_moment[2]) * hx.vec_mag(self.total_moment) / self.rho / (2.*self.rp)**5. / (self.w/2./np.pi)**2.
		self.CP = self.CL * 2. * np.pi
		print('Coefficient of thrust is {}'.format(self.CT))
		print('Coefficient of torque is {}'.format(self.CL))
		print('Coefficient of power is  {}\n\n'.format(self.CP))
	
	def pll(self, tol, max_it):
		print('\a\n\nAdvance ratio is {}\n\n'.format(self.J))
		self.renew_nu()
		if self.Ct_guess == -1.:
			pass
		else:
			self.first_guess()
		self.calc_gamma(tol, max_it)
		self.calc_forces_and_moments()
		self.calc_coef()
	
	def renew_nu(self):
		# calculate nu_ji
		self.nu = np.zeros((self.cp_total,self.cp_total, 3))
		print('Finished calculating velocity induced by horseshoe vortex #')
		for j in range(self.cp_total):
			a1 = (self.node[j,0,0] ** 2. + self.node[j,0,1] ** 2.) ** .5
			a2 = (self.node[j,1,0] ** 2. + self.node[j,1,1] ** 2.) ** .5
			phi = np.arctan2(self.cp[j,1],self.cp[j,0])
			#print('phi is {}'.format(phi))
			for i in range(self.cp_total):
				
				x = self.cp[i,0]
				y = self.cp[i,1]
				z = self.cp[i,2]
				r1 = np.array([x - self.node[j,0,0], y - self.node[j,1,0], 0.])
				r2 = np.array([x - self.node[j,0,1], y - self.node[j,1,1], 0.])
				
				
				#self.nu[i,j,:] = hx.nu(a1, a2, self.b, self.cp[i,:], self.n, self.m, self.p, self.rl, i!=j, phi)
				w = hf.calc_nu(self.nu[i,j,:],a1,a2,self.b[j,0],self.b[j,1],self.rl,self.cp[i,:],self.n,self.m,self.p,phi,'r',i==j)
				
			print(j+1, end=' ', flush=True)
	
	def renew_Vp(self):
		self.Vp = np.zeros((self.cp_total,3))
		for i in range(self.cp_total):
			self.Vp[i,:] = self.Vinf + np.cross(self.r[i,:],self.W)
	
	def local_iterative_pitch_solver(self, tol1, tol2, maxit):
		
		count = 1
		ea_max = 1.
		en = ea_max
		temp1 = self.Omega1
		
		while ea_max > tol1:
			
			eo = en
			print('Manual pitch set to {}'.format(self.manual_pitch))
			#print(self.b[:self.nb,:], '\n\n', self.b[self.nb:,:])
			print('\n\n\nPitch Iteration {}'.format(count))
			#input()
			self.pll(tol2, maxit)
			self.Ct_guess = -1.
			
			# individual pitch solvers per node
			temp = np.zeros([self.cp_total,2])
			for i in range(self.nb):
				for j in range(self.n_blades):
					k = i + j * self.nb
					
					if i == 0:
						
						temp[k,0] = self.V[k,2] / self.w * 2. * np.pi
						vavg = ( self.V[k,2] + self.V[k+1,2] ) / 2.
						temp[k,1] = vavg / self.w * 2. * np.pi
						
					elif i == self.nb - 1:
						
						temp[k,1] = self.V[k,2] / self.w * 2. * np.pi
						vavg = ( self.V[k,2] + self.V[k-1,2] ) / 2.
						temp[k,0] = vavg / self.w * 2. * np.pi
						
					else:
						
						vavg = (self.V[k-1,2] + self.V[k,2] ) / 2.
						temp[k,0] = vavg / self.w * 2. * np.pi
						vavg = (self.V[k,2] + self.V[k+1,2] ) / 2.
						temp[k,1] = vavg / self.w * 2. * np.pi
						
					
					#temp[k,0] = vavg / self.w * 2. * np.pi
					
					#if i < self.nb-1:
					#	temp[k,1] = self.V[k+1,2] / self.w * 2. * np.pi
					#else:
					#	temp[k,1] = temp[k,0]
			
			ea_max = 0.
			for i in range(self.cp_total):
				for j in range(2):
					ea = abs((temp[i,j] - self.b[i,j]) / temp[i,j])
					if ea > ea_max: ea_max = ea
			print('Max Approximate Error for the pitch solver is {}%\n\n'.format(ea_max*100.))
			count += 1
			en = ea_max
			if en >= eo:
				temp1 -= .1
			if temp1 < 0.:
				input('Local iterative pitch solver failed to converge!')
				break
			
			for i in range(self.cp_total):
				for j in range(2):
					self.b[i,j] += temp1 * (temp[i,j] - self.b[i,j])
					if self.b[i,j] < 0.: self.b[i,j] = .0
	
	def global_iterative_pitch_solver(self, tol1, tol2, maxit):
		
		count = 1
		ea = 1.
		
		while ea > tol1:
			
			print('Manual pitch set to {}'.format(self.manual_pitch))
			print('\n\n\nPitch Iteration {}'.format(count))
			print('Current global pitch is {}\n\n'.format(self.b[0,0]))
			#input()
			self.pll(tol2, maxit)
			self.Ct_guess = -1.
			
			# global pitch solvers per node
			avg_v = 0.
			for i in range(self.cp_total):
				avg_v += self.V[i,2]
			avg_v /= float(self.cp_total)
			temp = avg_v / self.w * 2. * np.pi
			
			ea = abs((temp - self.b[0,0]) / temp)
			print('Approximate Error for the global pitch solver is {}%\n\n'.format(ea*100.))
			count += 1
			
			self.b[:,:] += self.Omega1 * (temp - self.b[0,0])
			#if self.b[0,0] <= 0.: self.b[:,:] = .01
	
	
