import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri'#percorso dei file .txt
NOME_SPETTRO = 'Cs137_2.txt' #modificare con il nome del file
PATH = os.path.join(PATH, NOME_SPETTRO) 

counts = np.loadtxt(PATH, skiprows=12, max_rows=2048, unpack = True) #salta i commenti ed acquisice i conteggi dei canali 0-2047
channels = np.array([i for i in range(0, 2048)], dtype = float) #numero di canali

"""Convoluzione per trovare i picchi e le loro larghezze."""

convolved = np.zeros(len(counts))
convolved2 = np.zeros(len(counts))

"""Trova la delta ottimale per il filtro gaussiano di convoluzione, ovvero quella per cui 
la doppia convoluzione ha un minimo."""
delta_min = 1
delta_max = 40
min_convolv = np.array([])
for delta in range(delta_min, delta_max + 1):
    """Prima convoluzione con la derivata prima della gaussiana."""
    for j in range(len(counts)):
        convolved[j] = (-counts*((channels - j)/delta**2)*np.exp(-((channels - j)**2/(2*delta**2)))).sum()

    """Seconda convoluzione della convoluzione precedente con la derivata prima della gaussiana."""
    for j in range(len(counts)):
        convolved2[j] = (-convolved*((channels - j)/delta**2)*np.exp(-((channels - j)**2/(2*delta**2)))).sum()
    
    min_convolv = np.append(min_convolv, [np.amin(convolved2)])

delta = np.argmin(min_convolv) + 1
print(f'Delta ottimale = {delta}\n')

"""Prima convoluzione con la derivata prima della gaussiana."""
for j in range(len(counts)):
    convolved[j] = (-counts*((channels - j)/delta**2)*np.exp(-((channels - j)**2/(2*delta**2)))).sum()

"""Seconda convoluzione della convoluzione precedente con la derivata prima della gaussiana."""
for j in range(len(counts)):
    convolved2[j] = (-convolved*((channels - j)/delta**2)*np.exp(-((channels - j)**2/(2*delta**2)))).sum()


"""Trova tutti gli indici neg_ind per cui la doppia convoluzione è negativa e maggiore in modulo 
di una certa prominence (distanza dalla base), i picchi corrispondono
ai minimi della doppia convoluzione."""

[neg_ind, _] = find_peaks(-convolved2, prominence = 10)
print(f'Canali dei picchi: {neg_ind}\n')

"""Interpolazione lineare per trovare gli zeri di channels-convolved2."""
zc_indx = np.where(np.diff(np.sign(convolved2))) #indice prima di passare lo zero
conv_zero = np.array([]) #array che contiene gli indici degli zeri
for zc_i in zc_indx:
    t1 = channels[zc_i]
    t2 = channels[zc_i]
    a1 = convolved2[zc_i]
    a2 = convolved2[zc_i + 1]
    conv_zero = np.append((conv_zero), [(t1 + (0 - a1) * ((t2 - t1) / (a2 - a1)))])

print(f'Zeri della funzione di doppia convoluzione: {conv_zero[0:2*len(neg_ind)]}\n')

#conv_zero_left = conv_zero[:len(neg_ind) + 3:2]
#conv_zero_right = conv_zero[1:len(neg_ind) + 4:2]
#print(conv_zero_left, conv_zero_right)

sigma_left = np.array(np.sqrt(np.abs((conv_zero[:2*len(neg_ind):2] - neg_ind)**2 - 2*delta**2)))
sigma_right = np.array(np.sqrt(np.abs((conv_zero[1:2*len(neg_ind):2] - neg_ind)**2 - 2*delta**2)))

print(f'sigma_left: {sigma_left}\n')
print(f'sigma_right: {sigma_right}\n')

"""Metodo per trovare il background. m corrisponde alla FWHM del picco più largo, da trovare con
l'algoritmo per i picchi."""
media = 1600
sigma = 27
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
plt.title('Counts and background'+ ' ' + NOME_SPETTRO + ' ' + f'$\delta$ = {delta}')
plt.xlabel('Channels')
plt.ylabel('Counts')
#plt.plot(conv_zero, np.zeros(len(conv_zero)), 'o')
plt.plot(neg_ind, convolved2[neg_ind], marker = 'P', label = 'Picchi')
plt.plot(channels, convolved2, marker = 'o', label = 'Doppia convoluzione con $-\dfrac{x_{channel}-y}{\delta^2}e^{-\dfrac{(x_{channel}-y)^2}{2\delta^2}}$')
plt.plot(channels, counts, marker = 'o', color = 'b', label = 'Dati')
plt.plot(channels, k, marker = 'o', color = 'r', label = 'Fondo')
plt.plot(channels, net_counts, marker = 'o', color = 'skyblue', label = 'Net counts')
plt.minorticks_on()
plt.legend()
plt.show()
