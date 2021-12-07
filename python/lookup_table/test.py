import numpy as np
import lookup_table_reader as ltr
import helix_func as hf

lamb = float(input('lambda: '))
cyl = np.zeros(3)
cyl[0] = float(input('sigma: '))
cyl[1] = float(input('theta/pi: ')) * np.pi

v = ltr.read_table_quadratic(1.,lamb,cyl[0] * np.cos(cyl[1]), cyl[0] * np.sin(cyl[1]) )
x = ltr.read_table_linear(1.,lamb,cyl[0] * np.cos(cyl[1]), cyl[0] * np.sin(cyl[1]) )

w = hf.integrate_nondim(np.array([0.,0.,0.]), lamb, -1, cyl, 1000., 100., 1.2, 0., 'r')

#error = np.zeros(3)
error = abs((w - v) / w) * 100.
print(v,'lookup table quadratic')
print(w,'integrator')
print(error,'% error')
print()
error = abs((w - x) / w) * 100.
print(x,'lookup table linear')
print(w,'integrator')
print(error,'% error')
print()
