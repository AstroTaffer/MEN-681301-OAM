import numpy as np
from astropy.modeling import models
from scipy.optimize import curve_fit
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astroquery.xmatch import XMatch
from astroquery.vizier import Vizier
import pandas as pd

"""Columns
G = Vizier(catalog='I/355/gaiadr3', columns=['RA_ICRS', 'DE_ICRS', 'Gmag', 'e_Gmag'])
TM = Vizier(catalog='II/246/out', columns=['RAJ2000', 'DEJ2000', 'Kmag', 'Jmag', 'Hmag',
                                           'e_Kmag', 'e_Jmag', 'e_Hmag'])
MSX = Vizier(catalog='V/114/msx6_gp', columns=['RAJ2000', 'DEJ2000', 'B1', 'B2', 'A', 'C', 'D',
                                               'e_B1', 'e_B2', 'e_A', 'e_C', 'e_D'])
WISE = Vizier(catalog='II/311/wise', columns=['RAJ2000', 'DEJ2000', 'W1mag', 'W2mag', 'W3mag', 'W4mag',
                                              'e_W1mag', 'e_W2mag', 'e_W3mag', 'e_W4mag'])
PS = Vizier(catalog="II/349/ps1", columns=['RAJ2000', 'DEJ2000', 'gmag', 'rmag', 'imag', 'zmag', 'ymag',
                                           'e_gmag', 'e_rmag', 'e_imag', 'e_ymag'])
"""

def Planck(w, T, omega):
    b = models.BlackBody(temperature=T * u.K)
    F = b(w * u.nm).to(u.Jy / u.sr) * np.pi * u.sr * omega
    return F.value

coord = SkyCoord("2:49:8.40", "+60:34:12", unit=(u.hourangle, u.deg), frame='icrs')  # LP08
# coord = SkyCoord("11.798", "+85.244", unit=(u.deg, u.deg), frame='icrs')  # NGC188

G = Vizier(columns = ['RA_ICRS', 'DE_ICRS', 'Gmag', 'e_Gmag'])
G.ROW_LIMIT = -1
table_G = G.query_region(coord,
                        radius=10*u.arcmin,
                        catalog="I/355/gaiadr3")[0]

table_TM = XMatch.query(cat1=table_G,
                     cat2='vizier:II/246/out',
                     max_distance=0.5*u.arcsec,
                     colRA1='RA_ICRS', colDec1='DE_ICRS')

table_WISE = XMatch.query(cat1=table_TM,
                     cat2='vizier:II/311/wise',
                     max_distance=0.5*u.arcsec,
                     colRA1='RA_ICRS', colDec1='DE_ICRS')

table_PS1 = XMatch.query(cat1=table_WISE,
                     cat2='vizier:II/349/ps1',
                     max_distance=0.5*u.arcsec,
                     colRA1='RA_ICRS', colDec1='DE_ICRS')

ZP = np.array([3229, 1594, 1024, 666.7, 309.54, 171.79, 31.676, 8.3635, 3964.03, 3173.02, 2575.36, 2261.81, 2180.4])
wl = np.array([622, 1235, 1662, 2159, 3368, 4618, 12082, 22194, 481, 617, 752, 866, 962]) #cm

mag = [['Gmag', 'e_Gmag'], ['Jmag_2', 'e_Jmag_2'],['Hmag_2', 'e_Hmag_2'], ['Kmag_2', 'e_Kmag_2'],
       ['W1mag', 'e_W1mag'], ['W2mag', 'e_W2mag'], ['W3mag', 'e_W3mag'], ['W4mag', 'e_W4mag'],
       ['gmag', 'e_gmag'], ['rmag', 'e_rmag'], ['imag', 'e_imag'], ['zmag', 'e_zmag'], ['ymag', 'e_ymag']]

name = [['G_flux', 'e_G_flux', 'G_wl'],
         ['J_flux', 'e_J_flux', 'J_wl'],
         ['H_flux', 'e_H_flux', 'H_wl'],
         ['K_flux', 'e_K_flux', 'K_wl'],
         ['W1_flux', 'e_W1_flux', 'W1_wl'],
         ['W2_flux', 'e_W2_flux', 'W2_wl'],
         ['W3_flux', 'e_W3_flux', 'W3_wl'],
         ['W4_flux', 'e_W4_flux', 'W4_wl'],
         ['g_flux', 'e_g_flux', 'g_wl'],
         ['r_flux', 'e_r_flux', 'r_wl'],
         ['i_flux', 'e_i_flux', 'i_wl'],
         ['z_flux', 'e_z_flux', 'z_wl'],
         ['y_flux', 'e_y_flux', 'y_wl']]

t = pd.DataFrame()
t['RA_ICRS'] = table_PS1['RA_ICRS']
t['DE_ICRS'] = table_PS1['DE_ICRS']

for i in range(len(ZP)):
    for j in range(2):
        t[mag[i][j]] = table_PS1[mag[i][j]]
t.fillna(0., inplace=True)

for i in range(len(name)):
    t[name[i][0]] = ZP[i]*10**(-0.4*t[mag[i][0]]) # Jy
    t[name[i][1]] = t[mag[i][1]]*0.4*ZP[i]*10**(-0.4*t[mag[i][0]])*np.log(10) # Jy
    t[name[i][2]] = wl[i]*len(ZP) # cm

t = t.drop(columns=['Gmag', 'e_Gmag', 'Jmag_2', 'e_Jmag_2','Hmag_2', 'e_Hmag_2', 'Kmag_2', 'e_Kmag_2',
       'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag',
       'gmag', 'e_gmag', 'rmag', 'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag', 'ymag', 'e_ymag'])

print(t)

# Planck function
t = t.sort_values(by=['J_flux'], ascending=False)
t = t.reset_index(drop=True)

for i in range(10):
    flux = [t[name[j][0]][i] for j in range(len(wl))]
    flux_error = [t[name[j][1]][i] for j in range(len(wl))]
    popt, pcov = curve_fit(Planck, wl, flux, p0=(5000, 1e-20), bounds=([3000, 1e-23], [7000, 1e-18]))
    wl_data = np.linspace(450, 23e3, 500)
    fig, ax = plt.subplots(layout='constrained')
    ax.errorbar(wl[0], flux[0], yerr=flux_error[0], fmt="o", markersize=5, c='red', label='GaiaDR3')
    ax.errorbar([wl[1], wl[2], wl[3]], [flux[1], flux[2], flux[3]], yerr=[flux_error[1], flux_error[2], flux_error[3]], fmt="o", markersize=5, c='green', label='2MASS')
    ax.errorbar([wl[4], wl[5], wl[6], wl[7]], [flux[4], flux[5], flux[6], flux[7]],
                yerr=[flux_error[4], flux_error[5], flux_error[6], flux_error[7]], fmt="o", markersize=5, c='orange', label='WISE')
    ax.errorbar([wl[8], wl[9], wl[10], wl[11], wl[12]], [flux[8], flux[9], flux[10], flux[11], flux[12]],
                yerr=[flux_error[8], flux_error[9], flux_error[10], flux_error[11], flux_error[12]], fmt="o", markersize=5, c='yellow', label='Pan-Starrs')
    ax.plot(wl_data, Planck(wl_data, *popt))
    ax.set_xlabel(r'lg{$\lambda$}, nm')
    ax.set_ylabel('lgF, Jy')
    ax.loglog()
    ax.legend()
    plt.show()