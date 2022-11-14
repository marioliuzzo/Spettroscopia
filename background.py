import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri'#percorso dei file .txt
NOME_SPETTRO = 'Na22_1.txt' #modificare con il nome del file
PATH = os.path.join(PATH, NOME_SPETTRO) 

counts = np.loadtxt(PATH, skiprows=12, max_rows=2048, unpack = True) #salta i commenti ed acquisice i conteggi dei canali 0-2047
channels = np.array([i for i in range(0, 2048)], dtype = float) #numero di canali

#for i in range(media - 3*sigma, media + 3*sigma + 1):
#    AREA_NETTA = net_counts[i].sum()

#print(f'{AREA_NETTA:.3f}')

"""Convoluzione per trovare i picchi e le loro larghezze."""

delta = 23 #deviazione standard del filtro gaussiano
convolved = np.zeros(len(counts))
convolved2 = np.zeros(len(counts))

"""Prima convoluzione con la derivata prima della gaussiana."""
for j in range(len(counts)):
    convolved[j] = (-counts*((channels - j)/delta**2)*np.exp(-((channels - j)**2/(2*delta**2)))).sum()

"""Seconda convoluzione della convoluzione precedente con la derivata prima della gaussiana."""
for j in range(len(counts)):
    convolved2[j] = (-convolved*((channels - j)/delta**2)*np.exp(-((channels - j)**2/(2*delta**2)))).sum()

"""Trova tutti gli indici neg_ind per cui la doppia convoluzione è negativa e maggiore in modulo 
di una certa prominence (distanza dalla base) i picchi corrispondono
ai minimi della doppia convoluzione."""
[neg_ind, _] = find_peaks(-convolved2, prominence = 20)
print(neg_ind)





"""Metodo per trovare il background. m corrisponde alla FWHM del picco più largo, da trovare con
l'algoritmo per i picchi."""
media = 1600
sigma = 26
m = int(np.floor(2.35*sigma))
k = np.array([i for i in counts])
z = np.zeros(len(counts))
for p in range(1, m+1):
    for i in range(p-1, len(counts)-p):
        z[i] = min(k[i], 0.5*(k[i-p] + k[i+p]))
    for i in range(p-1, len(counts)-p):
        k[i] = z[i]

net_counts = np.array([i for i in counts]) - k

NOME_SPETTRO = NOME_SPETTRO.replace('_2.txt', '')
NOME_SPETTRO = NOME_SPETTRO.replace('_1.txt', '')
plt.title('Counts and background'+ ' ' + NOME_SPETTRO)
plt.xlabel('Channels')
plt.ylabel('Counts')
plt.plot(neg_ind, convolved2[neg_ind], marker = 'P', label = 'Picchi')
plt.plot(channels, convolved2, marker = 'o', linestyle = ' ', label = 'Doppia convoluzione con $-\dfrac{x_{channel}-y}{\delta^2}e^{-\dfrac{(x_{channel}-y)^2}{2\delta^2}}$')
plt.plot(channels, counts, marker = 'o', color = 'b', label = 'Dati')
#plt.plot(channels, k, marker = 'o', color = 'r', label = 'Fondo')
#plt.plot(channels, net_counts, marker = 'o', label = 'Net counts')
plt.minorticks_on()
plt.legend()
plt.show()
