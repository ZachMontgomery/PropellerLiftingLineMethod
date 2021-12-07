import numpy as np
import helix_func as hf
import helix as hx
import time

a1 = .1
a2 = .2
b = 1.
n = 200.
m = 20.
p = 1.2
rl = -1.
i = 1
j = 10
x = 0.
y = 1.
z = 0.
phi = 10. * np.pi / 180.
control = np.array([x,y,z])


py = np.zeros(3)
start_p = time.time()
py = hx.nu(a1, a2, b, control, n, m, p, rl, i!=j, phi)
finish_p = time.time()


fo = np.zeros(3)
start_f = time.time()
fo = hf.calc_nu(fo,a1,a2,b,rl,control,n,m,p,phi,'r',i==j)
finish_f = time.time()

time_p = finish_p - start_p
time_f = finish_f - start_f
dif = time_p - time_f



print('Python results are  ',py)
print('Fortran results are ',fo,'\n\n')

print('Python  ran for {:8.4f} seconds'.format(time_p))
print('Fortran ran for {:8.4f} seconds'.format(time_f))
print('Difference is   {:8.4f} seconds'.format(dif))
print('Fortran is {:12.4f}% faster'.format(dif / time_f * 100.))
