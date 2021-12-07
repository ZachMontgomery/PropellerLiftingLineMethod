import numpy as np
import matplotlib.pyplot as plt
import re

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




bet = open('blade_loading.txt','r')
pll = open('pll_blade_loading.txt','r')

bet_res = 25
pll_res = 25

zeta_bet = np.zeros(bet_res)
zeta_pll = np.zeros(pll_res)

y_bet = np.zeros(bet_res)
y_pll = np.zeros(pll_res)

z_bet = np.zeros(bet_res)
z_pll = np.zeros(pll_res)

for i in range(bet_res):
	temp = myfunc(bet)
	zeta_bet[i] = float(temp[0])
	y_bet[i] = float(temp[1])
	z_bet[i] = float(temp[2])
	

for i in range(pll_res):
	temp = myfunc(pll)
	zeta_pll[i] = float(temp[0])
	y_pll[i] = float(temp[2])
	z_pll[i] = float(temp[3])
	

plt.figure()
plt.plot(zeta_bet,y_bet,'--k')
plt.plot(zeta_pll,y_pll,'k')
plt.title('Comparison')
plt.xlabel('Blade Radius')
plt.ylabel('Y Force lbf (fighting rotation)')
plt.axis([0.,1.,0.,.065])
plt.legend(['Blade Element Theory','Lifting-Line'])
plt.grid()

plt.figure()
plt.plot(zeta_bet,-z_bet,'--k')
plt.plot(zeta_pll,-z_pll,'k')
plt.title('Comparison')
plt.xlabel('Blade Radius')
plt.ylabel('Z Force lbf (thrust)')
plt.axis([0.,1.,0.,.2])
plt.legend(['Blade Element Theory','Lifting-Line'])
plt.grid()

plt.show()
