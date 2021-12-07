import blade_mod as bm
import numpy as np
import json

g = json.load(open('propeller.json'))



res = g['propeller']['blades']['control points']

J = g['scenario']['advance ratio']
z = np.zeros(res)
y = np.zeros(res)
f = open('bet_blade_loading.txt', 'w')


zeta, y, z, CT, Cl, twist, CT_p = bm.blade_element_theory(J,res)
for i in range(res):
	line = '{0}		{1}		{2}		{3}\n'.format(zeta[i], y[i], z[i], twist[i])
	f.write(line)

f.write('{}		{}'.format(CT, Cl))

f.close()

print('coefficient of thrust from BET is {}'.format(CT))
print('coefficient of thrust from P   is {}'.format(CT_p))

import matplotlib.pyplot as plt
plt.figure()
plt.plot(zeta,y)
plt.plot(zeta,-z)
plt.legend(['Drag','Thrust'])
plt.show()
