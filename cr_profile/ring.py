from ring_function import *

from astromodels import *

import astropy.units as u

#RA and Dec of OB2 association

starRA = 308.3  
starDec = 41.3167 

cocoon_index = -2.65

rings = []
spectra = []
sources = []



rings.append( Disk_on_sphere() )

for i in range(0,4):
  rings.append( Ring_on_sphere() )

for i in range(0,5):
  spectra.append( Powerlaw() )
  sources.append( ExtendedSource("Ring%d" %i,  spectral_shape = spectra[i], spatial_shape=rings[i]  ) )
  
for i, ring in enumerate( rings ):
  inner = 0.6*(i)
  outer = 0.6*(1+i)
  ring.lon0=starRA*u.degree
  ring.lat0=starDec*u.degree
  ring.lon0.fix=True
  ring.lat0.fix=True
  
  if i>0:
    ring.radius_inner = inner*u.degree
    ring.radius_inner.fix = True
    ring.radius_inner.bounds = (0.2,inner)*u.degree
    ring.radius_outer = outer*u.degree
    ring.radius_outer.fix = True
    ring.radius_outer.bounds = (0.2,outer)*u.degree

  else:
    ring.radius = outer * u.degree
    ring.radius.fix = True
    ring.radius.bounds = (0.2,outer)*u.degree
   

for spectrum in spectra:
  spectrum.piv = 4.2*u.TeV
  spectrum.piv.fix = True 
  spectrum.K = 3.5e-14 / u.cm**2 / u.s / u.TeV
  spectrum.K.bounds = (1e-19, 1e-8) / u.cm**2 / u.s / u.TeV
  spectrum.K.fix = False
  #spectrum.K.set_uninformative_prior(Log_uniform_prior)   
  spectrum.index = cocoon_index
  spectrum.index.fix = False
  spectrum.index.bounds = (-4, 0)

  model = Model(*sources ) 
  
model.save("results/test_likelihood_0.6_free.yml")
