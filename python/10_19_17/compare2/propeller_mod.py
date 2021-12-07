import numpy as np
import json
import helix as hx
import helix_func as hf
import blade_shape_pll as bs

class propeller:
	
	def __init__(self, filename):
		
		f = json.load(open(filename))
		
		# number of blades
		self.n_blades = f['propeller']['blades']['number']
		# number of control points per blade
		self.nb = f['propeller']['blades']['control points']
		# total number of control points
		self.cp_total = self.n_blades * self.nb
		# propeller radius
		self.rp = f['propeller']['tip radius']
		# hub radius
		self.rh = f['propeller']['hub radius']
		# right-handed or left-handed propeller
		if f['propeller']['direction of rotation'] == 'right-handed':
			self.rl = -1.
		elif f['propeller']['direction of rotation'] == 'left-handed':
			self.rl = 1.
		else:
			input('\n\n{1}\n{2}Error in {0}!\n{1}\n\nDirection of rotation should be \'right-handed\' or \'left-handed\''.format(filename,100*'*',int((90-len(filename))/2)*' '))
		# angular velocity of propeller
		self.w = f['propeller']['RPMs'] / 60. * 2. * np.pi
		# number of loops in helices
		self.n = f['numerics']['helix']['loops']
		# downstream distance
		self.zh = f['numerics']['helix']['downstream distance']
		# avg number of steps per loop
		self.m = f['numerics']['helix']['steps']
		# power clustering factor
		self.p = f['numerics']['helix']['power factor']
		# relaxation factor for newtons method
		self.Omega = f['numerics']['relaxation factor']['newtons method']
		# relaxation factor for pitch iteration
		self.Omega1 = f['numerics']['relaxation factor']['pitch solver']
		# fluid density
		self.rho = f['scenario']['density']
		# chord-line pitch length
		self.Kc = f['propeller']['pitch'] / 2. / self.rp
		# manual pitch
		self.manual_pitch = f['propeller']['manual pitch'] * np.pi / 180.
		# contraction variables
		self.c1 = f['numerics']['helix']['contraction distance'] / 4.
		self.c2 = f['numerics']['helix']['contraction percentage']
		# freestream velocity vs advance ratio
		if 'advance ratio' in f['scenario']:
			self.J = f['scenario']['advance ratio']
			self.vinf = self.J * self.w * self.rp / np.pi
		elif 'freestream velocity' in f['scenario']:
			self.vinf = f['scenario']['freestream velocity']
			self.J = np.pi * self.vinf / self.w / self.rp
		else:
			input('\n\n{1}\n{2}Error in {0}!\n{1}\n\n{0} must contain a value for either \'freestream velocity\' or \'advance ratio\'.\n\nAdvance ratio value overrides freestream velocity value if both are given!'.format(filename,100*'*',int((90-len(filename))/2)*' '))
		# pitch
		if 'Ct guess' in f['scenario']:
			self.update_pitch(Ct = f['scenario']['Ct guess'])
		else:
			self.update_pitch()
		# freestream velocity vector
		self.Vinf = np.array([0.,0.,self.vinf])
		# angular velocity vector
		self.W = np.array([0.,0.,self.rl*self.w])
		
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
			self.twist[i] = bs.twist(self.cp[i,0], self.manual_pitch, self.Kc)
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
		
		# calculate dl
		self.dl = np.zeros((self.cp_total,3))
		for i in range(self.cp_total):
			self.dl[i,:] = self.rl*(self.node[i,0,:] - self.node[i,1,:])
		
		# calculate dA
		self.dA = np.zeros((self.cp_total))
		for i in range(self.cp_total):
			self.dA[i] = (self.c[i,0] + self.c[i,1]) / 2. * hx.vec_mag(self.dl[i,:]) #(hx.vec_mag(self.node[i,1,:]) - hx.vec_mag(self.node[i,0,:]))
	
	def update_pitch(self, **arg):
		
		# Guess for CT
		if 'Ct' in arg:
			self.Ct_guess = arg['Ct']
		else:
			self.Ct_guess = .05 #float(input('\n\nAdvance ratio is {}, Enter guess for CT: '.format(self.J)))
		# Vg
		if self.vinf == 0.:
			self.Vg = 2. * self.w * self.rp * np.sqrt(self.Ct_guess / np.pi ** 3.) * 0.8
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
		
		#~ first guess without small angle approximations
		b = np.zeros(self.cp_total)
		A = np.zeros((self.cp_total,self.cp_total))
		for i in range(self.cp_total):
			b[i] = self.Vg
			for j in range(self.cp_total):
				A[i,j] = self.nu[i,j,2]
		
		# #~ first guess using small angle approximations
		# #~ right side matrix
		# b = np.zeros((self.cp_total))
		# for i in range(self.cp_total):
			# b[i] = hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.Vp[i,:],self.un[i,:]) * self.dA[i]
		# #~ calculate A matrix
		# A = np.zeros((self.cp_total, self.cp_total))
		# for i in range(self.cp_total):
			# for j in range(self.cp_total):
				# if i == j:
					# A[i,j] = 2. * hx.vec_mag(np.cross(self.Vp[i,:],self.dl[i,:])) - hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.nu[i,j,:],self.un[i,:]) * self.dA[i]
				# else:
					# A[i,j] = - hx.vec_mag(self.Vp[i,:]) ** 2. * 2. * np.pi * np.dot(self.nu[i,j,:],self.un[i,:]) * self.dA[i]
		
		self.gamma = np.linalg.solve(A,b)
		# print(self.gamma)
		# input()
	
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
				temp -= .05
				print('Relaxed Newton\'s Method further to promote convergence. Relaxation factor is {}'.format(temp))
				if temp <= 1.e-10:
					input()
					break
			
			print('Iteration = {:<3d}  Max residual = {:6.1e}'.format(i,e))
			
			if i >= max_it:
				print('Newtons method failed to converge in {} iterations!'.format(i))
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
			
			spanwise_vel = np.dot(self.V[i,:],self.dl[i,:]) / hx.vec_mag(self.dl[i,:]) ** 2. * self.dl[i,:]
			V_temp = self.V[i,:] - spanwise_vel
			
			
			self.dF[i,:] = self.rho * self.gamma[i] * np.cross(self.V[i,:],self.dl[i,:]) + .5 * self.rho * hx.vec_mag(V_temp) * (self.c[i,0] + self.c[i,1]) / 2. * bs.Cd(self.aoa[i]) * hx.vec_mag(self.dl[i,:]) * V_temp
			self.total_force += self.dF[i,:]
			
			self.us[i,:] = np.cross(self.ua[i,:],self.un[i,:])
			self.Cm[i] = -.1
			self.dM[i,:] = -1./8.*self.rho*hx.vec_mag(self.V[i,:])**2.*self.Cm[i]*(self.c[i,0]+self.c[i,1])**2.*hx.vec_mag(self.dl[i,:])*self.us[i,:]
			self.total_moment += np.cross(self.r[i,:],self.dF[i,:]) + self.dM[i,:]
		print('\n\nThe total force is  {}'.format(self.total_force))
		print('The total moment is {}\n\n'.format(self.total_moment))
	
	def calc_coef(self):
		self.CT = -self.total_force[2] / self.rho / (2.*self.rp)**4. / (self.w/2./np.pi)**2.
		self.CL = -self.rl * self.total_moment[2] / self.rho / (2.*self.rp)**5. / (self.w/2./np.pi)**2.
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
		# ## calculate nu_ji
		# self.nu = np.zeros((self.cp_total,self.cp_total, 3))
		# print('Finished calculating velocity induced by horseshoe vortex #')
		# for j in range(self.cp_total):
			# a1 = (self.node[j,0,0] ** 2. + self.node[j,0,1] ** 2.) ** .5
			# a2 = (self.node[j,1,0] ** 2. + self.node[j,1,1] ** 2.) ** .5
			# phi = np.arctan2(self.cp[j,1],self.cp[j,0])
			# #print('phi is {}'.format(phi))
			# for i in range(self.cp_total):
				
				# #~ x = self.cp[i,0]
				# #~ y = self.cp[i,1]
				# #~ z = self.cp[i,2]
				# #~ r1 = np.array([x - self.node[j,0,0], y - self.node[j,1,0], 0.])
				# #~ r2 = np.array([x - self.node[j,0,1], y - self.node[j,1,1], 0.])
				
				# if self.b[j,0] <= 0.:
					# self.b[j,0] = .01
				# if self.b[j,1] <= 0.:
					# self.b[j,1] = .01
				
				# self.n0 = self.zh / self.b[j,0]
				# self.n1 = self.zh / self.b[j,1]
				
				# w = hf.calc_nu(self.nu[i,j,:],a1,a2,self.b[j,0],self.b[j,1],self.rl,self.cp[i,:],self.n0,self.n1,self.m,self.p,phi,'r',i==j)
				
			# print(j+1, end=' ', flush=True)
		# print('\n\n')
		
		# initialize nu_ji
		self.nu = np.zeros((self.cp_total,self.cp_total, 3))
		print('Finished calculating velocity induced by horseshoe vortex #')
		# calculate nu on entire first blade for every control point
		for j in range(self.nb):
			a1 = (self.node[j,0,0] ** 2. + self.node[j,0,1] ** 2.) ** .5
			a2 = (self.node[j,1,0] ** 2. + self.node[j,1,1] ** 2.) ** .5
			phi = np.arctan2(self.cp[j,1],self.cp[j,0])
			for i in range(self.cp_total):
				if self.zh > 0.:
					if self.b[j,0] <= 0.:
						self.b[j,0] = .1
					if self.b[j,1] <= 0.:
						self.b[j,1] = .1
					self.n0 = self.zh / self.b[j,0]
					self.n1 = self.zh / self.b[j,1]
				else:
					self.n0 = self.n
					self.n1 = self.n
				w = hf.calc_nu(self.nu[i,j,:],a1,a2,self.b[j,0],self.b[j,1],self.rl,self.cp[i,:],self.n0,self.n1,self.m,self.p,phi,'r',i==j)
			print(j+1, end=' ', flush=True)
		ang_spacing = 2. * np.pi / float(self.n_blades)
		# loop through for each additional blade
		for k in range(1,self.n_blades):
			# loop through the horseshoes on the blade
			for jj in range(self.nb):
				# translate jj to the nu matrix
				j = jj + k * self.nb
				# loop through all the control points
				for ii in range(self.cp_total):
					# translate i to the new blade
					i = ii + k * self.nb
					if i >= self.cp_total: i -= self.cp_total
					# calculate rotation
					theta = float(k) * ang_spacing
					x,y,z = self.nu[ii,jj,:]
					xx = x * np.cos(theta) - y * np.sin(theta)
					yy = y * np.cos(theta) + x * np.sin(theta)
					# set nu value
					self.nu[i,j,:] = [xx,yy,z]
				print(j+1, end=' ', flush=True)
		print('\n\n')
		
		# for j in range(self.cp_total):
			# for i in range(self.cp_total):
				# er = hx.vec_mag(self.nu[i,j,:] - nu[i,j,:])
				# if er >= 1.e-13:
					# print('{:6.1e} error at control point {} and horseshoe {}'.format(er,i,j)) #self.nu[i,j,:],'\n',nu[i,j,:],
		# input()
	
	def renew_Vp(self):
		self.Vp = np.zeros((self.cp_total,3))
		for i in range(self.cp_total):
			self.Vp[i,:] = self.Vinf + np.cross(self.r[i,:],self.W)
	
	def local_iterative_pitch_solver(self, tol1, tol2, maxit):
		
		def interpolate(self, y1, y2, y, x1, x2):
			return (y - y1) * (x2 - x1) / (y2 - y1) + x1
		
		count = 1
		ea_max = 1.
		en = ea_max
		temp1 = self.Omega1
		
		while ea_max > tol1:
			
			eo = en
			tol = 1.e-10
			for i in range(self.nb):
				for j in range(2):
					for k in range(self.n_blades):
						if abs(self.b[i,j] - self.b[i+k*self.nb,j]) > tol:
							print(self.b[i,j])
							print(self.b[i+k*self.nb,j])
							print(i)
							print(i+k*self.nb)
							print(abs(self.b[i,j] - self.b[i+k*self.nb,j]),'\n')
							self.b[i+k*self.nb,j] = self.b[i,j]
			
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
						
						vel = interpolate(self,self.cp[i,0],self.cp[i+1,0],self.node[i,0,0],self.V[i,2],self.V[i+1,2])
						temp[k,0] = vel / self.w * 2. * np.pi
						vel = interpolate(self,self.cp[i,0],self.cp[i+1,0],self.node[i,1,0],self.V[i,2],self.V[i+1,2])
						temp[k,1] = vel / self.w * 2. * np.pi
						
					elif i == self.nb - 1:
						
						vel = interpolate(self,self.cp[i-1,0],self.cp[i,0],self.node[i,1,0],self.V[i-1,2],self.V[i,2])
						temp[k,1] = vel / self.w * 2. * np.pi
						vel = interpolate(self,self.cp[i-1,0],self.cp[i,0],self.node[i,0,0],self.V[i-1,2],self.V[i,2])
						temp[k,0] = vel / self.w * 2. * np.pi
						
					else:
						
						vel = interpolate(self,self.cp[i-1,0],self.cp[i,0],self.node[i,0,0],self.V[i-1,2],self.V[i,2])
						temp[k,0] = vel / self.w * 2. * np.pi
						vel = interpolate(self,self.cp[i,0],self.cp[i+1,0],self.node[i,1,0],self.V[i,2],self.V[i+1,2])
						temp[k,1] = vel / self.w * 2. * np.pi
						
					
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
			print('Max Approximate Error for the local pitch solver is {:6.1e}\n\n'.format(ea_max))
			count += 1
			en = ea_max
			# if en >= eo:
				# temp1 -= .001
			if temp1 < 0.:
				input('Local iterative pitch solver failed to converge!')
				break
			
			for i in range(self.cp_total):
				for j in range(2):
					self.b[i,j] += temp1 * (temp[i,j] - self.b[i,j])
					# if self.b[i,j] <= 0.:
						# self.b[i,j] = .01
						# input('manually set pitch')
	
	def global_iterative_pitch_solver(self, tol1, tol2, maxit):
		
		count = 1
		ea = 1.
		
		while ea > tol1:
			
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
			print('Approximate Error for the global pitch solver is {:6.1e}\n\n'.format(ea))
			count += 1
			
			self.b[:,:] += self.Omega1 * (temp - self.b[0,0])
	
	def streamline_solver(self):
		
		import matplotlib.pyplot as plt
		
		k = int(self.nb / 2)
		oo = np.pi / 2.
		rr = self.cp[k,0]
		zz = 0.
		count = 10
		
		plt.ion()
		f = open('streamline_data.csv', 'w')
		
		while zz <= 20.:
			res = 10
			dzeta = 1.e-4
			
			r = np.zeros(res)
			o = np.zeros(res)
			z = np.zeros(res)
			y = np.zeros(res)
			x = np.zeros(res)
			vz = np.zeros(res-1)
			vx = np.zeros(res-1)
			vy = np.zeros(res-1)
			vr = np.zeros(res-1)
			vo = np.zeros(res-1)
			
			r[0] = rr
			o[0] = oo
			z[0] = zz
			
			for k in range(res):
				x[k] = r[k] * np.cos(o[k])
				y[k] = r[k] * np.sin(o[k])
				
				nu = np.zeros((self.cp_total,3))
				V = np.array([0.,0.,self.vinf])
				if k < res-1:
					for j in range(self.cp_total):
						a1 = (self.node[j,0,0] ** 2. + self.node[j,0,1] ** 2.) ** .5
						a2 = (self.node[j,1,0] ** 2. + self.node[j,1,1] ** 2.) ** .5
						phi = np.arctan2(self.cp[j,1],self.cp[j,0])
						#~ #nu[j,:] = ch.nu(a1,a2,self.b[j,0],self.b[j,1],self.n,self.c1,self.c2,phi,self.rl,x[k],y[k],z[k],self.m,self.p,self.cp_total,j)
						w = hf.calc_nu(nu[j,:],a1,a2,self.b[j,0],self.b[j,1],self.rl,[x[k],y[k],z[k]],self.n,self.m,self.p,phi,'r',False)
						V[:] += self.gamma[j] * nu[j,:]
						#print(j, end = ' ', flush = True)
					vz[k] = V[2]
					
					vr[k] = V[0] * np.cos(o[k]) + V[1] * np.sin(o[k])
					vo[k] = V[1] * np.cos(o[k]) - V[0] * np.sin(o[k]) - self.rl * r[k] * self.w
					
					vx[k] = vr[k] * np.cos(o[k]) - vo[k] * np.sin(o[k])
					vy[k] = vr[k] * np.sin(o[k]) + vo[k] * np.cos(o[k])
					
					count += 1
					
					if count >= 10:
						f.write('{},{},{},{},{},{},{},{},{},{}\n'.format(x[k],y[k],z[k],r[k],o[k],vx[k],vy[k],vz[k],vr[k],vo[k]))
						#print(z[k])
						count = 0
					
					r[k+1] = r[k] + dzeta * vr[k]
					o[k+1] = o[k] + dzeta * vo[k]
					z[k+1] = z[k] + dzeta * vz[k]
				
				
				print('{:5.1f}%'.format(float(k+1)/float(res)*100.))
			
			#~ vx[res-1] = vx[res-2]
			#~ vy[res-1] = vy[res-2]
			#~ vz[res-1] = vz[res-2]
			
			rr = r[res-1]
			oo = o[res-1]
			zz = z[res-1]
			
			plt.figure(1)
			plt.plot(x,y,'k')
			plt.xlabel('x')
			plt.ylabel('y')
			plt.figure(2)
			plt.plot(z,x,'b')
			plt.plot(z,y,'g')
			plt.xlabel('z')
			plt.ylabel('x or y')
			plt.figure(3)
			plt.plot(z,r,'k')
			plt.xlabel('z')
			plt.ylabel('r')
			plt.figure(4)
			plt.plot(z[:res-1],vx,'b')
			plt.plot(z[:res-1],vy,'g')
			plt.xlabel('z')
			plt.ylabel('vx or vy')
			plt.figure(5)
			plt.plot(z[:res-1],vz,'k')
			plt.xlabel('z')
			plt.ylabel('vz')
			plt.draw()
			plt.pause(.1)
			#~ #input('press enter to continue')
		plt.show()
		input()
		
	
	
