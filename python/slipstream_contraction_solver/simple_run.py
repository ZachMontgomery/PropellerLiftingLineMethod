import time
import numpy as np
from propeller_mod import *

start = time.time()

prop = propeller('propeller.json')

prop.global_iterative_pitch_solver(1.e-1, 1.e-12, 100)

close = time.time()
print('\n\nProgram ran for {:8.3f} seconds'.format(close - start))


#print('  no freestream velocity {}'.format(prop.b[0,0] / 4. / prop.rp * np.sqrt( np.pi / prop.CT ) ))
#print('with freestream velocity {}'.format((prop.b[0,0]**2.*prop.w*np.pi-prop.b[0,0]*prop.vinf*np.pi**2.)/16./prop.CT/prop.w/prop.rp**2.))


prop.iterate_slipstream_contraction_variables()
