import numpy as np
from propeller_mod import *
import blade_mod as bm
import matplotlib.pyplot as plt


f = open('initial_guess_test.txt','w')


prop = propeller('propeller.json')

tol_pitch = 1.e-5
tol_nm = 1.e-12
maxit = 100

prop.global_iterative_pitch_solver(tol_pitch, tol_nm, maxit)


f.close()

