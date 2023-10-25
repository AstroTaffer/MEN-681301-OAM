import numpy as np
from astropy.io import ascii
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
from astroquery.xmatch import XMatch
from astroquery.vizier import Vizier


coord = SkyCoord("2:49:8.40", "+60:34:12", unit=(u.hourangle, u.deg), frame='icrs')  # LP08
# coord = SkyCoord("11.798", "+85.244", unit=(u.deg, u.deg), frame='icrs')  # NGC188
w = u.Quantity(10, u.arcmin)
h = u.Quantity(10, u.arcmin)

"""table = XMatch.query(cat1='vizier:I/355/gaiadr3',
                     cat2='vizier:II/246/out',
                     max_distance=5 * u.arcsec, colRA1='ra',
                     colDec1='dec')
print(table)"""

result_1 = Vizier.query_region(coord,
                        width=w, height=h,
                        catalog=["I/355/gaiadr3", "II/246/out", "V/114/msx6_gp", "II/311/wise", "II/349/ps1"])
# result = Vizier.query_object("NGC188", catalog=["I/355/gaiadr3", "II/246/out", "V/114/msx6_gp", "II/311/wise", "II/349/ps1"])
print(result_1)

G = Vizier(catalog='I/355/gaiadr3', columns=['RA_ICRS', 'DE_ICRS', 'Gmag'])
TM = Vizier(catalog='II/246/out', columns=['RAJ2000', 'DEJ2000', 'Kmag', 'Jmag', 'Hmag'])
MSX = Vizier(catalog='V/114/msx6_gp', columns=['RAJ2000', 'DEJ2000', 'B1', 'B2', 'A', 'C', 'D'])
WISE = Vizier(catalog='II/311/wise', columns=['RAJ2000', 'DEJ2000', 'W1mag', 'W2mag', 'W3mag', 'W4mag'])
PS = Vizier(catalog="II/349/ps1", columns=['RAJ2000', 'DEJ2000', 'gmag', 'rmag', 'imag', 'zmag', 'ymag'])

result_G = G.query_region(coord, radius=10*u.arcmin)
result_TM = TM.query_region(coord, radius=10*u.arcmin)
result_WISE = WISE.query_region(coord, radius=10*u.arcmin)
result_PS = PS.query_region(coord, radius=10*u.arcmin)

table = result_G[0]
RA_GAIA = table['RA_ICRS']
DE_GAIA = table['DE_ICRS']

