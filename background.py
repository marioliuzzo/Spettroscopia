import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri'#percorso dei file .txt
NOME_SPETTRO = 'Cs137_2.txt' #modificare con il nome del file
PATH = os.path.join(PATH, NOME_SPETTRO) 

M = 10 #numero canali prima e dopo il picco preso dal fitgauss
counts = np.loadtxt(PATH, skiprows=12, max_rows=2048, unpack = True) #salta i commenti ed acquisice i conteggi dei canali 0-2047
channels = np.array([i for i in range(0, 2048)], dtype = float) #numero di canali

channels1 = np.array([channels[i] for i in range(790, 940)], dtype = float) #canali vicino al picco, da n a n_max-1
counts1 = np.array([counts[i] for i in range(790, 940)], dtype = float)


counts2 = np.array([counts[i] for i in range(790-M, 791)], dtype = float)#canali prima del picco
counts3 = np.array([counts[i] for i in range(939, 939 + M + 1)], dtype = float)#canali dopo il picco

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
SIGMA = np.sqrt(C + Fondo*len(channels1)/2)
print(f'Area netta = {AREA_NETTA:.3f} +- {SIGMA:.3f}')
'''
CA = counts[820]
CB = counts[900]
Fondo = (CA + CB)*0.5*len(channels1)
AREA_TOT = counts1.sum()
AREA_NETTA = AREA_TOT - Fondo
print(AREA_NETTA)