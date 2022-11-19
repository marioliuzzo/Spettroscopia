import fitgauss as fit
import numpy as np
"""
------------------------------------------------------------
Metodo per trovare il fondo con i canali adiacenti al picco
------------------------------------------------------------
"""

channels = fit.channels
counts = fit.counts

channels1 = np.array([channels[i] for i in range(70, 110)], dtype=float)
counts1 = np.array([counts[i] for i in range(70, 110)], dtype=float)

init_values = init_values = [94., 5., 500000., 2000.]

F = fit.FitGauss(channels1, counts1, init_values)
iterator = iter(F)
elem = np.array([])
while True:
    try:
        # Get next element from TeamIterator object using iterator object
        elem = np.append(elem, [next(iterator)])
    except StopIteration:
        break
mu0 = elem[0]
sigma0 = elem[1]
A0 = elem[2]
B0 = elem[3]
dm, dsigma, dA, dB = np.sqrt(F.covm.diagonal())

#canali a 3 sigma dal picco, trovato con il fit
A = int(np.floor(mu0 - 3*sigma0))
B = int(np.floor(mu0 + 3*sigma0))

AREA_TOT = 0
for i in range(A, B + 1):
    AREA_TOT += counts[A]

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

N = channels[B] - channels[A]
FONDO = 0.5*(mu_i + mu_f)*N
VAR_FONDO = FONDO * N/(2*m)
SIGMA_FONDO = np.sqrt(VAR_FONDO)

AREA_NETTA = AREA_TOT - FONDO
SIGMA_NETTA = np.sqrt(VAR_AREA_TOT + VAR_FONDO)

FWHM = 2.35*sigma0
RISOLUZIONE = FWHM/mu0
print(f'Risoluzione = {RISOLUZIONE:.3f}\n')

print(f'fondo = {FONDO:.3f} +- {SIGMA_FONDO:.3f}\n')
print(f'Area netta = {AREA_NETTA:.3f} +- {SIGMA_NETTA:.3f}\n')
print(f'media = {mu0:.3f} +- {dm:.3f}\n')
print(f'sigma = {sigma0:.3f} +- {dsigma:.3f}\n')
print(f'A = {A0:.3f} +- {dA:.3f}\n')
print(f'B = {B0:.3f} +- {dB:.3f}\n')
