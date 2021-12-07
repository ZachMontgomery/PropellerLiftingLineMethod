import numpy as np
import re
import json

numeric_const_pattern = r"""
     [-+]? # optional sign
     (?:
         (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
         |
         (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
     )
     # followed by optional exponent part if desired
     (?: [Ee] [+-]? \d+ ) ?
     """
rx = re.compile(numeric_const_pattern, re.VERBOSE)
del numeric_const_pattern

def myfunc(f):
	dummy = rx.findall(f.readline())
	return dummy

pll = open('pll_blade_loading.txt', 'r')
bet = open('bet_blade_loading.txt', 'r')

f = json.load(open('propeller.json'))

res = f['propeller']['blades']['control points']
blades = f['propeller']['blades']['number']
J = f['scenario']['advance ratio']
tol = 1e14

zeta_pll = np.zeros(res)
zeta_bet = np.zeros(res)
thrust_pll = np.zeros(res)
thrust_bet = np.zeros(res)
drag_pll = np.zeros(res)
drag_bet = np.zeros(res)
twist_pll = np.zeros(res)
twist_bet = np.zeros(res)



for i in range(res):
	zeta_pll[i], drag_pll[i], thrust_pll[i], twist_pll[i] = myfunc(pll)
	zeta_bet[i], drag_bet[i], thrust_bet[i], twist_bet[i] = myfunc(bet)
	
	if abs(zeta_pll[i] - zeta_bet[i]) > tol:
		print('Error in zeta is {} at control point {}'.format(abs(zeta_pll[i] - zeta_bet[i]),i))
	elif abs(twist_pll[i] - twist_bet[i]) > tol:
		print('Error in twist is {} at control point {}'.format(abs(twist_pll[i] - twist_bet[i]),i))

CT_pll, CL_pll = myfunc(pll)
CT_bet, CL_bet = myfunc(bet)

CT_pll = float(CT_pll)
CT_bet = float(CT_bet)
CL_pll = float(CL_pll)
CL_bet = float(CL_bet)

import matplotlib.pyplot as plt

f, axarr = plt.subplots(3, sharex=True, figsize = (10,30))


axarr[2].plot(zeta_pll,-thrust_pll)
axarr[2].plot(zeta_bet,-thrust_bet)
axarr[2].legend(['PLL','BET'])
axarr[2].set_ylabel('Thrust')
axarr[2].set_xlabel('Blade Span')
axarr[2].set_title('PLL CT = {:10.8f}\nBET CT = {:10.8f}\n{:4.1f}% Difference'.format(CT_pll,CT_bet,abs((CT_pll - CT_bet) / CT_pll) * 100.))

axarr[1].plot(zeta_pll,drag_pll)
axarr[1].plot(zeta_bet,drag_bet)
axarr[1].legend(['PLL','BET'])
axarr[1].set_ylabel('Drag')
axarr[1].set_title('PLL CL = {:10.8f}\nBET CL = {:10.8f}\n{:4.1f}% Difference'.format(CL_pll,CL_bet,abs((CL_pll-CL_bet)/CL_pll)*100.))

axarr[0].plot(zeta_pll,twist_pll*180./np.pi)
axarr[0].plot(zeta_bet,twist_bet*180./np.pi)
axarr[0].legend(['PLL','BET'])
axarr[0].set_ylabel('Twist (deg)')
axarr[0].set_title('Advance Ratio = {:4.1f}   blades = {}'.format(J,blades))



# plt.figure()
# plt.plot(zeta_pll,-thrust_pll)
# plt.plot(zeta_bet,-thrust_bet)
# plt.legend(['PLL','BET'])
# plt.ylabel('Thrust')
# plt.xlabel('Blade Span')
# plt.title('PLL CT = {}\nBET CT = {}'.format(CT_pll,CT_bet))

# plt.figure()
# plt.plot(zeta_pll,drag_pll)
# plt.plot(zeta_bet,drag_bet)
# plt.legend(['PLL','BET'])
# plt.ylabel('Drag')
# plt.xlabel('Blade Span')
# plt.title('PLL CL = {}\nBET CL = {}'.format(CL_pll,CL_bet))

# plt.figure()
# plt.plot(zeta_pll,twist_pll*180./np.pi)
# plt.plot(zeta_bet,twist_bet*180./np.pi)
# plt.legend(['PLL','BET'])
# plt.ylabel('Twist Angle (deg)')
# plt.xlabel('Blade Span')

plt.show()
