import helix as hx
import numpy as np

#v = np.zeros(3)
b = .5
n = 100.
m = 50.
p = 1.2
rl = 1.
phi = 0.

#v = hx.rk4_v(a, b, x, y, z, n, m, p, rl, phi)
#w = hx.rk4_v(a, b, x, -y, z, n, m, p, -rl, phi)

#print(v,'right-handed')
#print(w,'left-handed')

xp = 0.
yp = 0.
zp = 0.
j1 = .1
j2 = .2

r1 = np.array([xp-j1,yp,zp])
r2 = np.array([xp-j2,yp,zp])

nu = hx.rk4_v(j1, b, xp, yp, zp, n, m, p, rl, phi)
print(nu, 'h1')
nu = hx.rk4_v(j2, b, xp, yp, zp, n, m, p, rl, phi)
print(nu, 'h2')

nu = hx.nu(j1, j2, b, xp, yp, zp, n, m, p, rl, r1, r2, 1, 2, phi)

print(nu, 'nu')
