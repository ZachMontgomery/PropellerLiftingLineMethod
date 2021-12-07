import numpy as np
from propeller_mod import *
import matplotlib.pyplot as plt
import json
import blade_mod as bm

filename = 'propeller.json'
f = open('C_vs_J_bet.txt','w')


res_grid = 5
nb = 5

for k in range(res_grid):
	
	nb *= 2
	
	g = json.load(open(filename))
	g['propeller']['blades']['control points'] = nb
	new_file = json.dumps(g, indent = 4)
	ff = open(filename,'w')
	ff.write(new_file)
	ff.close()
	
	# prop = propeller(filename)
	
	# tol_pitch = 1.e-5
	# tol_nm = 1.e-12
	# maxit = 300
	
	
	
	
	res_j = 11
	
	
	for i in range(res_j):
		
		J = 1. - float(i) / float(res_j-1)
		# J = 0.
		
		# prop.update_advance_ratio(J=J, Ct = 0.0)
		
		# prop.global_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
		# prop.pll(tol_nm, maxit)
		
		CT, CP, CL = bm.blade_element_theory(J,nb)
		
		# line = '{0}		{1}		{2}		{3}		{4}\n'.format(J, prop.CT, prop.CP, prop.CL, prop.nb)
		line = '{0}		{1}		{2}		{3}		{4}\n'.format(J, CT, CP, CL, nb)
		f.write(line)
	
	
	


f.close()

