import numpy as np
from astropy.io import ascii
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt

# coord = SkyCoord("2:49:8.40", "+60:34:12", unit=(u.hourangle, u.deg), frame='icrs')  # LP08
coord = SkyCoord("11.798", "+85.244", unit=(u.deg, u.deg), frame='icrs')  # NGC188
w = u.Quantity(20, u.arcmin)
h = u.Quantity(20, u.arcmin)
Gaia.ROW_LIMIT = 500
r = Gaia.query_object_async(coordinate=coord, width=w, height=h)

fig1, ax1 = plt.subplots(dpi=100)
ax1.set_xlabel(r'$\alpha$, deg')
ax1.set_ylabel(r'$\delta$, deg')
ax1.scatter(r["ra"], r["dec"], c=r["phot_g_mean_mag"], s=15)
plt.show()

fig2, ax2 = plt.subplots(dpi=100)
ax2.set_xlabel(r'$\alpha$, deg')
ax2.set_ylabel(r'$\delta$, deg')
plt.quiver(r["ra"].filled(0), r["dec"].filled(0), r["pmra"].filled(0), r["pmdec"].filled(0), r["phot_g_mean_mag"])
plt.show()
