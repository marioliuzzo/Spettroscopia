import logging
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

"""
======================================================================================================
Fit gaussiano nell'intorno dei picchi, per ricavare media (posizione del canale), deviazione standard,
ampiezza ed offset.
======================================================================================================
"""

PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri'  # percorso dei file .txt
NOME_SPETTRO = 'Am241_1.txt'  # modificare con il nome del file
PATH = os.path.join(PATH, NOME_SPETTRO)

# mette i risultati del fit nel file NOME_SPETTROlog.txt
logging.basicConfig(filename=NOME_SPETTRO.replace('txt', '_bckg_log.txt'),
                    level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

# salta i commenti ed acquisice i conteggi dei canali 0-2047
counts = np.loadtxt(PATH, skiprows=12, max_rows=2048, unpack=True)
channels = np.array([i for i in range(0, 2048)],
                    dtype=float)  # numero di canali

# canali vicino al picco, da n a n_max-1
channels1 = np.array([channels[i] for i in range(443, 523)], dtype=float)
counts1 = np.array([counts[i] for i in range(443, 523)], dtype=float)


def gaussiana(x, mu, sigma, A, B):
    """Funzione per fit gaussiano channels-counts. A è l'ampiezza della gaussiana
    B è l'offset."""
    return A*(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mu)/sigma)**2) + B


def log_results(channel, mu, dmu, sigma, dsigma, A, dA, B, dB):
    """Inserisce i risultati del fit nel file NOME_SPETTRO_bckg.log.txt ."""
    logging.info(f'Range di canali: {channel[0]}-{channel[-1]}\n')
    logging.info(f'media = {mu:.3f} +- {dmu:.3f}\n')
    logging.info(f'dev. std = {sigma:.3f} +- {dsigma:.3f}\n')
    logging.info(f'A = {A:.3f} +- {dA:.3f}\n')
    logging.info(f'B = {B:.3f} +- {dB:.3f}\n')


def plot_results(channels, counts, mu, sigma, A, B, NOME_SPETTRO):
    """Plotta i risultati del fit nell'intervallo counts1, channels1."""
    plt.plot(channels, gaussiana(channels, mu, sigma, A, B), color='red')
    plt.plot(channels, counts, marker='o')
    NOME_SPETTRO = NOME_SPETTRO.replace('_1.txt', '')
    plt.title('Channels vs counts' + ' ' + NOME_SPETTRO)
    plt.xlabel('Channels')
    plt.ylabel('Counts')
    plt.minorticks_on()
    plt.show()


class Fit_iterator:
    """Iteratore di FitGauss."""

    def __init__(self, classe):
        self.classe = classe
        self.i = 0

    def __next__(self):
        if self.i < len(self.classe.Fit()[0]):
            result = self.classe.Fit()[0][self.i]
            self.i += 1
            return result
        raise StopIteration


class FitGauss:
    """Classe per il fit gaussiano."""

    def __init__(self, x, y, init):
        self.init = init
        self.x = x
        self. y = y
        self.pars = np.array([])
        self.covm = np.array([], [])

    def Fit(self):
        self.pars, self.covm = curve_fit(gaussiana, self.x, self.y, self.init)
        return self.pars, self.covm

    def __iter__(self):
        return Fit_iterator(self)


if __name__ == '__main__':

    init_values = [482., 21., 500000., 2000.]
    F = FitGauss(channels1, counts1, init_values)
    iterator = iter(F)
    elem = np.array([])
    while True:
        try:
            # Get next element from TeamIterator object using iterator object
            elem = np.append(elem, [next(iterator)])
        except StopIteration:
            break
    mu0 = elem[3]
    sigma0 = elem[2]
    A0 = elem[1]
    B0 = elem[0]
    dm, dsigma, dA, dB = np.sqrt(F.covm.diagonal())
    print(f'{mu0:.3f} + {dm:.3f}')
#    log_results(channels1, mu0, dmu, sigma0, dsigma, A0, dA, B0, dB)
#    plot_results(channels1, counts1, mu0, sigma0, A0, B0, NOME_SPETTRO)
