import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""
================================================================================================
Fit lineare per la calibrazione canali-energia, con le sorgenti di Am241, Ba133, Cs137_2, Co61, Na22
================================================================================================
"""


Energy, Channel, dChannel = np.loadtxt('linear.txt', unpack = True)
dEnergy = np.zeros(len(Energy))

def linear_fit(x, a, b):
    return a*x + b

init_values = [1., 10.]
pars, covm = curve_fit(linear_fit, Energy, Channel, init_values, dChannel)
a0, b0 = pars
da, db = np.sqrt(covm.diagonal())

chisq = (((Channel - linear_fit(Energy, *pars))/dChannel)**2).sum()
ndof = len(Energy) - 2

print(f'a = {a0:.3f} +- {da:.3f}')
print(f'b = {b0:.3f} +- {db:.3f}')
print(f'chisq/ndof = {chisq:.3f}/{ndof}')
print(covm)

plt.errorbar(Energy, Channel, dChannel, dEnergy, marker = 'o', linestyle = '')
x = np.linspace(Energy[0], Energy[-1], 100)
plt.plot(x, linear_fit(x, *pars), color = 'r')
plt.xlabel('Energy [KeV]')
plt.ylabel('Channels')
plt.show()
