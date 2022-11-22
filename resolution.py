import numpy as np
import linearfit as linear
import adj_chan as adj

"""
Calcolo della risoluzione energetica dalla media e deviazione standard date
dal fit gaussiano allo spettro al netto di fondo e continuum, in adj_chan.
La media trovata viene convertita in energia [keV] attraverso la calibrazione,
data dal fit in linearfit.
"""

FWHM = adj.FWHM
sigma = adj.sigma1
dsigma = adj.dsigma1

pars = linear.pars
covm = linear.covm
da, db = np.sqrt(covm.diagonal())
correlation = covm[0][1]/(da*db)  # correlazione fra a e b

# canali che si vogliono convertire in energia, da usare per il calcolo della risoluzione
# è il mu1 di adj_chan, media dello spettro al netto di fondo e continuum
x = adj.mu1

ener_calibration = pars[0]*x + pars[1]  # energia ottenuta dalla calibrazione
dener_calibration = np.sqrt(
    covm[0][0]*x**2 + covm[1][1] + 2*x*correlation*da*db)  # incertezza sull'energia
# considerando la correlazione fra i parametri
RISOLUZIONE = FWHM/ener_calibration
DRISOLUZIONE = 2.35*np.sqrt((dsigma/ener_calibration) **
                            2 + (sigma/(ener_calibration**2))**2*dener_calibration**2)

print(f'\n Risoluzione energetica = {100*RISOLUZIONE:.3f} +- {100*DRISOLUZIONE:.3f} %\n')
