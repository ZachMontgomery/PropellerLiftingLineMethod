import time
import numpy as np
from propeller_mod import *

start = time.time()

prop = propeller('propeller.json')

res = 31
J = np.zeros(res)
CT = np.zeros(res)

for i in range(res):
	J[i] = float(i) / float(res-1) * 0.2
	prop.update_advance_ratio(J=J[i], Ct=0.0)
	prop.pll(1.e-12, 100)
	CT[i] = prop.CT

close = time.time()
print('\n\nProgram ran for {:8.3f} seconds'.format(close - start))

import matplotlib.pyplot as plt
plt.figure()
plt.plot(J,CT)
plt.show()

