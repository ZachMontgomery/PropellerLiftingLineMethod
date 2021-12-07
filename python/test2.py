from propeller_mod import *
import numpy as np

prop = propeller('input.txt')

prop.renew_Vp()

print(prop.b, prop.Vinf)

print('Advance ratio is {}'.format(prop.J))
input()

prop.renew_nu()

prop.pll(1.e-8, 5000)

print('CT is ',prop.CT)
print('CL is ',prop.CL)
print('CP is ',prop.CP)




