import numpy as np
from propeller_mod import *
import blade_mod as bm
import matplotlib.pyplot as plt


f = open('pll_C_vs_J.txt','w')
g = open('bet_C_vs_J.txt','w')

res = 51

prop = propeller('propeller.json')

tol_pitch = 1.e-5
tol_nm = 1.e-12
maxit = 100







for i in range(res):
	
	J = 1. - float(i) / float(res-1)
	
	prop.update_advance_ratio(J=J)
	
	prop.global_iterative_pitch_solver(tol_pitch, tol_nm, maxit)
	
	CT, CP, CL = bm.blade_element_theory(J,prop.nb)[3:6]
	
	line = '{0}		{1}		{2}		{3}\n'.format(J, prop.CT, prop.CP, prop.CL)
	f.write(line)
	
	line = '{0}		{1}		{2}		{3}\n'.format(J, CT, CP, CL)
	g.write(line)
	


f.close()
g.close()

