import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy_healpix import HEALPix
import astropy.units as u
from scipy.interpolate import griddata

def eq2gal(ra, dec):
    c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    l, b = c.galactic.l.deg, c.galactic.b.deg
    return l, b

F = np.load('I_355_tgextmap.npz')
lst = F.files

# ra = 201+26/60+47/3600
# dec = -(47+28/60+51/3600)

# ra = 11.798
# dec = 85.244
ra = 254.77817888
dec = -42.70233451
ra = ra*u.deg
dec = dec*u.deg

"""l, b = eq2gal(ra, dec)
l = l*u.deg
b = b*u.deg"""
"""for item in lst:
    print(item)
    print(F[item])"""
level = 9
HPId = F['HPId']
A0 = F['A0']
m = F['HPlevel'] == level

# hp = HEALPix(nside=2**level, order='nested', frame='galactic')
hp = HEALPix(nside=2**level, order='nested', frame='icrs')
radius = 15 * u.deg
# coords = SkyCoord(l=l, b=b, frame='galactic')
coords = SkyCoord(ra, dec, frame='icrs')
search = hp.cone_search_skycoord(coords, radius=radius)
ext = A0[np.isin(HPId, search)]
common = HPId[np.isin(HPId, search)]

cd = hp.healpix_to_skycoord(common)
"""l1 = np.unique(cd.l.value)
b1 = np.unique(cd.b.value)
X, Y = np.meshgrid(l1, b1)"""
ra1 = np.unique(cd.ra.value)
dec1 = np.unique(cd.dec.value)
X, Y = np.meshgrid(ra1, dec1)

beta = griddata(np.transpose(np.array([cd.ra.value, cd.dec.value])), np.array(ext), (X, Y), method='nearest')
im = plt.imshow(beta, cmap=plt.cm.viridis, aspect='auto', extent=[X.min(), X.max(), Y.min(), Y.max()],
                origin='lower')
plt.colorbar(im)
plt.show()




