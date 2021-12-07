import helix_func as hx
import numpy as np


nu = np.zeros(3)
a1 = .1
a2 = .2
b = .5
rl = 1.
control = np.array([0., 0., 0.])
n = 100.
m = 50.
p = 1.2
phi = 0.
int_type = 'r'

nu = hx.integrate(nu, a1, b, rl, control, n, m, p, phi, int_type)
#print(nu)
nu = hx.integrate(nu, a2, b, rl, control, n, m, p, phi, int_type)
#print(nu)

nu = hx.calc_nu(nu, a1, a2, b, rl, control, n, m, p, phi, int_type, False)
print(nu)

