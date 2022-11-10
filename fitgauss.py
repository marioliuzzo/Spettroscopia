import os 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



counts = np.loadtxt('Cs137_1.txt', skiprows=12, max_rows=2048, unpack = True)
channels = np.array([i for i in range(0, 2048)], dtype = float)

def gaussiana(x, mu, sigma, A, B):
    """Funzione per fit gaussiano channels-counts."""
    return A*(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mu)/sigma)**2) + B

channels1 = np.array([channels[i] for i in range(449, 511)], dtype = float)
counts1 = np.array([counts[i] for i in range(449, 511)], dtype = float)

init_values = [500., 21., 3000., 500.]
pars, covm = curve_fit(gaussiana, channels1, counts1, init_values)

mu0, sigma0, A0, B0 = pars
dmu, dsigma, dA, dB = np.sqrt(covm.diagonal())
#dy = np.ones(len(counts1))
#chisq = (((counts1-gaussiana(channels1, mu0, sigma0, A0, B0))/dy)**2).sum()
#ndof = len(channels1)-4
print(f'media = {mu0} +- {dmu}')
print(f'dev. std = {sigma0} +- {dsigma}')
print(f'A = {A0} +- {B0}')
print(f'B = {B0} +- {dB}')
#print(f'chisq/ndof = {chisq}/{ndof}')


plt.plot(channels1, counts1, marker = 'o')
plt.plot(channels1, gaussiana(channels1, mu0, sigma0, A0, B0), color = 'red')
plt.title('Channels vs counts')
plt.xlabel('Channels')
plt.ylabel('Counts')
plt.show()