import os
import matplotlib.pyplot as plt
import fitgauss as fit
import numpy as np
"""
------------------------------------------------------------
Metodo per trovare il fondo con i canali adiacenti al picco.
------------------------------------------------------------
"""

NOME_BCKG = 'fondo_54437.txt'
PATH_BCKG = os.path.join('C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri/', NOME_BCKG)

background = np.loadtxt(PATH_BCKG, skiprows=12, max_rows=2048, unpack=True)

NOME_SPETTRO = fit.NOME_SPETTRO.replace('_2.txt', '')
NOME_SPETTRO = fit.NOME_SPETTRO.replace('_1.txt', '')

#live time di acquisizione del background e dei radionuclidi
LIVE_TIME_BCKG = 54437
LIVE_TIME = {'Am241': 271, 'Ba133': 233,
             'Co60': 344, 'Na22': 1517, 'Cs137': 194}

#il background è riscalato con il rapporto fra il live time della
#misura dello spettro ed il live time della misura del fondo senza
#sorgente
background = background * LIVE_TIME[NOME_SPETTRO]/LIVE_TIME_BCKG

channels = fit.channels
counts = fit.counts

channels1 = fit.channels1
counts1 = fit.counts1

init_values = fit.init_values

F = fit.FitGauss(channels1, counts1, init_values)
risultati = fit.risultati(F)
mu0 = risultati[0]
sigma0 = risultati[1]
A0 = risultati[2]
B0 = risultati[3]
dm, dsigma, dA, dB = np.sqrt(F.covm.diagonal())

#canali a 3 sigma dal picco, trovato con il fit
A = int(np.floor(mu0 - 3*sigma0))
B = int(np.floor(mu0 + 3*sigma0))

A = 445
AREA_TOT = 0
for i in range(A, B):
    AREA_TOT += counts[i]

VAR_AREA_TOT = AREA_TOT
SIGMA_AREA_TOT = np.sqrt(VAR_AREA_TOT)

print(f'Area totale = {AREA_TOT:.3f} +- {SIGMA_AREA_TOT:.3f}\n')

#numero di canali da mediare prima di A e dopo B
m = 5

#media dei conteggi degli m canali prima di A e dopo B
mu_i = 0
mu_f = 0
for i in range(m+1):
    mu_i += counts[A-i]
    mu_f += counts[B+i]
mu_i = mu_i/m
mu_f = mu_f/m

N = int(channels[B] - channels[A])
AREA_CONTINUUM = 0.5*(mu_i + mu_f)*N
VAR_FONDO = AREA_CONTINUUM * N/(2*m)
SIGMA_FONDO = np.sqrt(VAR_FONDO)

for i in range(A, B):
    BCKG = background[i].sum()

#area ottenuta sottraendo continuum e fondo ai dati
AREA_NETTA = AREA_TOT - AREA_CONTINUUM - BCKG
SIGMA_NETTA = np.sqrt(VAR_AREA_TOT + VAR_FONDO + BCKG)

def retta(x, y, A, B):
    """Interpolazione lineare fra i punti A e B del continuum."""
    return ((x - A)/(B - A))*(y[B] - y[A]) + y[A]

x = np.linspace(int(channels[A]), int(channels[B]), N)
retta = retta(x, counts, A, B)

channels_restrict = np.array([])#canali ristretti fra A e B, per plottare la retta
counts_restrict = np.array([])#conteggi ristretti fra A e B, per il fit allo spettro netto
spettro_netto = np.array([])#spettro a cui si sottraggono fondo e continuum


for i in range(A, B):
    spettro_netto = np.append(spettro_netto, counts[i] - retta[i - A] - background[i])

for i in range(A, B):
    channels_restrict = np.append(channels_restrict, channels[i]) 
    counts_restrict = np.append(counts_restrict, counts[i])

init_values = [mu0, sigma0, A0, B0]
F1 = fit.FitGauss(channels_restrict, counts_restrict, init_values)
risultati = fit.risultati(F1)
mu1 = risultati[0]
sigma1 = risultati[1]
A1 = risultati[2]
B1 = risultati[3]

dm1, dsigma1, dA1, dB1 = np.sqrt(F1.covm.diagonal())
FWHM = 2.35*sigma1

print(f'Area continuum = {AREA_CONTINUUM:.3f} +- {SIGMA_FONDO:.3f}\n')
print(f'Area netta = {AREA_NETTA:.3f} +- {SIGMA_NETTA:.3f}\n')
#print(f'media = {mu0:.3f} +- {dm:.3f}\n')
#print(f'sigma = {sigma0:.3f} +- {dsigma:.3f}\n')
#print(f'A = {A0:.3f} +- {dA:.3f}\n')
#print(f'B = {B0:.3f} +- {dB:.3f}\n')
print(f'media netta = {mu1:.3f} +- {dm1:.3f}\n')
print(f'sigma netta = {sigma1:.3f} +- {dsigma1:.3f}\n')

plt.xlabel('Channels [UA]')
plt.ylabel('Counts [UA]')
plt.plot(channels, counts, marker='o', color='b', label='Dati')
plt.plot(channels, background, marker = 'o', label = 'Fondo true')
plt.plot(x, retta, marker = 'o', label = 'Retta fra i punti a $\pm$ $3\sigma$ dal picco')
plt.plot(channels_restrict, spettro_netto, marker = 'o', label = 'Spettro al netto di fondo e continuum')
plt.minorticks_on()
plt.legend()
plt.show()
