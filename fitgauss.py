import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

if __name__ == '__main__':
    PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri'#percorso dei file .txt
    NOME_SPETTRO = 'Am241_1.txt' #modificare con il nome del file
    PATH = os.path.join(PATH, NOME_SPETTRO) 

    #mette i risultati del fit nel file NOME_SPETTROlog.txt
    logging.basicConfig(filename = NOME_SPETTRO.replace('txt', '_bckg_log.txt'), level = logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

    counts = np.loadtxt(PATH, skiprows=12, max_rows=2048, unpack = True) #salta i commenti ed acquisice i conteggi dei canali 0-2047
    channels = np.array([i for i in range(0, 2048)], dtype = float) #numero di canali

    channels1 = np.array([channels[i] for i in range(69, 110)], dtype = float) #canali vicino al picco, da n a n_max-1
    counts1 = np.array([counts[i] for i in range(69, 110)], dtype = float)

    def gaussiana(x, mu, sigma, A, B):
        """Funzione per fit gaussiano channels-counts."""
        return A*(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mu)/sigma)**2) + B

    init_values = [91., 6., 500000., 2000.]
    pars, covm = curve_fit(gaussiana, channels1, counts1, init_values)

    mu0, sigma0, A0, B0 = pars
    dmu, dsigma, dA, dB = np.sqrt(covm.diagonal())

    logging.info(f'Range di canali: {channels1[0]}-{channels1[-1]}\n')  
    logging.info(f'media = {mu0:.3f} +- {dmu:.3f}\n')
    logging.info(f'dev. std = {sigma0:.3f} +- {dsigma:.3f}\n')
    logging.info(f'A = {A0:.3f} +- {B0:.3f}\n')
    logging.info(f'B = {B0:.3f} +- {dB:.3f}\n')


    plt.plot(channels1, gaussiana(channels1, mu0, sigma0, A0, B0), color = 'red')
    plt.plot(channels1, counts1, marker = 'o')
    NOME_SPETTRO = NOME_SPETTRO.replace('_1.txt', '')
    plt.title('Channels vs counts' + ' ' + NOME_SPETTRO)
    plt.xlabel('Channels')
    plt.ylabel('Counts')
    plt.show()
