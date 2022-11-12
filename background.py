import os
import logging
import numpy as np
import matplotlib.pyplot as plt

PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri'#percorso dei file .txt
NOME_SPETTRO = 'Cs137_2.txt' #modificare con il nome del file
PATH = os.path.join(PATH, NOME_SPETTRO) 

M = 10 #numero canali prima e dopo il picco preso dal fitgauss
counts = np.loadtxt(PATH, skiprows=12, max_rows=2048, unpack = True) #salta i commenti ed acquisice i conteggi dei canali 0-2047
channels = np.array([i for i in range(0, 2048)], dtype = float) #numero di canali

channels1 = np.array([channels[i] for i in range(786, 943)], dtype = float) #canali vicino al picco, da n a n_max-1
#counts1 = np.array([counts[i] for i in range(786, 943)], dtype = float)


counts2 = np.array([counts[i] for i in range(786-M, 787)], dtype = float)#canali prima del picco
counts3 = np.array([counts[i] for i in range(942, 942 + M + 1)], dtype = float)#canali dopo il picco

'''
Mu_i = 1/M * (counts2.sum())
Mu_f = 1/M * (counts3.sum())
Fondo = 0.5 * (Mu_i + Mu_f) * len(channels1)
print(Fondo)
C = counts1.sum()
print(C)
print(Mu_i)
print(Mu_f)
AREA_NETTA = C - Fondo
SIGMA = np.sqrt(C + Fondo*len(channels1)*0.5)
print(f'Area netta = {AREA_NETTA:.3f} +- {SIGMA:.3f}')
'''
media = 91
sigma = 26
m = int(np.floor(2.35*sigma))
k = np.array([i for i in counts])
z = np.zeros(len(counts))
for p in range(1, m+1):
    for i in range(p-1, len(counts)-p):
        z[i] = min(k[i], 0.5*(k[i-p] + k[i+p]))
    for i in range(p, len(counts)-p):
        k[i] = z[i]

net_counts = np.array([i for i in counts]) - k

for i in range(media - 3*sigma, media + 3*sigma):
    AREA_NETTA = net_counts[i].sum()

print(f'{AREA_NETTA}')

NOME_SPETTRO = NOME_SPETTRO.replace('_2.txt', '')
NOME_SPETTRO = NOME_SPETTRO.replace('_1.txt', '')
plt.title('Counts and background'+ ' ' + NOME_SPETTRO)
plt.xlabel('Channels')
plt.ylabel('Counts')
plt.plot(channels, counts, marker = 'o', color = 'b', label = 'Dati')
plt.plot(channels, k, marker = 'o', color = 'r', label = 'Fondo')
plt.plot(channels, net_counts, marker = 'o', label = 'Net counts')
plt.minorticks_on()
plt.legend()
plt.show()
