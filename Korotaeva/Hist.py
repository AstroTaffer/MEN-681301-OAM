import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import ascii
from astroquery.vizier import Vizier
from astroquery.skyview import SkyView

Gaia = Vizier(columns=["pmRA", "pmDE", "Plx", "Gmag"])
Gaia.ROW_LIMIT = -1
Gaia = Gaia.query_region("NGC 5139", radius=20 * u.arcmin, catalog="I/355/gaiadr3")[0]

pmRA = Gaia['pmRA'][np.where((Gaia["pmRA"] > -20) & (Gaia["pmRA"] < 10))]
pmDE = Gaia['pmDE'][np.where((Gaia["pmDE"] > -20) & (Gaia["pmDE"] < 10))]
Plx = Gaia['Plx'][np.where((Gaia["Plx"] > -10) & (Gaia["Plx"] < 10))]
Gmag = Gaia['Gmag']

fig, ax = plt.subplots(2, 2, figsize=(10, 10))

ax[0][0].hist(pmRA, bins=30, color='#56B4E9', edgecolor='black')
ax[0][0].set_xlabel(r'$\mu_\alpha$, mas/yr', size=14)
ax[0][0].set_ylabel('N', size=14)

ax[0][1].hist(pmDE, bins=30, color='#56B4E9', edgecolor='black')
ax[0][1].set_xlabel(r'$\mu_\delta$, mas/yr', size=14)

ax[1][0].hist(Plx, bins=30, color='#56B4E9', edgecolor='black')
ax[1][0].set_xlabel(r'$\pi$, mas', size=14)
ax[1][0].set_ylabel('N', size=14)

ax[1][1].hist(Gmag, bins=30, color='#56B4E9', edgecolor='black')
ax[1][1].set_xlabel('G, mag', size=14)
plt.show()