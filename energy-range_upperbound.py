from hawc_hal import HAL, HealpixConeROI, HealpixMapROI
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
mpl.rc("font", family="serif", size=14)

import warnings

from threeML.plugin_prototype import PluginPrototype
from astromodels.core.sky_direction import SkyDirection
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# This disable momentarily the printing of warnings, which you might get
# if you don't have the Fermi ST or other software installed
with warnings.catch_warnings():

    warnings.simplefilter("ignore")

    from threeML import *

#Load the data and detector response
maptree = "maptree.hd5"
response = "response.hd5"

#Spectral description
spectrum = Cutoff_powerlaw()
shape    = Gaussian_on_sphere()

source   = ExtendedSource("veritasellipse",
                          spatial_shape=shape,
                          spectral_shape=spectrum)

fluxUnit = 1. / (u.TeV * u.cm**2 * u.s)

# Set spectral parameters (do it after the source definition to make sure
# the units are handled correctly)
spectrum.K = 1.34e-13 * fluxUnit
spectrum.K.fix = True

spectrum.piv = 4.9 * u.TeV
spectrum.piv.fix = True


spectrum.index = -1.9
spectrum.index.fix = True

spectrum.xc = 33 * u.TeV
spectrum.xc.fix = True



# Set spatial parameters (do it after the source definition to make sure
# the units are handled correctly)
shape.lon0 = 307.89 * u.degree
shape.lon0.fix = True

shape.lat0 = 41.58 * u.degree
shape.lat0.fix = True

shape.sigma = 0.27 * u.degree
shape.sigma.fix = True
shape.sigma.max_value = 0.27 * u.degree


spectrum3 = Powerlaw()* StepFunction()
shape3    = Gaussian_on_sphere()

source3   = ExtendedSource("HAWCcocoon",
                          spatial_shape=shape3,
                          spectral_shape=spectrum3)

fluxUnit = 1. / (u.TeV * u.cm**2 * u.s)

# Set spectral parameters (do it after the source definition to make sure
# the units are handled correctly)
spectrum3.K_1 = 3.5e-13 * fluxUnit
spectrum3.K_1.fix = False
spectrum3.K_1.bounds = (1e-16 * fluxUnit, 1e-10 * fluxUnit)
#spectrum3.K.set_uninformative_prior(Log_uniform_prior)    

spectrum3.piv_1 = 4.2 * u.TeV
spectrum3.piv_1.fix = True


spectrum3.index_1 = -2.20
spectrum3.index_1.fix = False
spectrum3.index_1.bounds = (-4, 0)
#spectrum3.index.set_uninformative_prior(Uniform_prior)


spectrum3.lower_bound_2 = 0.1 * u.TeV
spectrum3.lower_bound_2.fix = True

spectrum3.upper_bound_2 = 100 * u.TeV
spectrum3.upper_bound_2.fix = False
spectrum3.upper_bound_2.bounds = (1 * u.TeV, 500. * u.TeV)

spectrum3.value_2 = 1.
spectrum3.value_2.fix = True

# Set spatial parameters (do it after the source definition to make sure
# the units are handled correctly)
shape3.lon0 = 307.65 * u.degree
shape3.lon0.fix = True

shape3.lat0 = 40.93 * u.degree
shape3.lat0.fix = True

shape3.sigma = 2.13 * u.degree
shape3.sigma.fix = True
shape3.sigma.max_value = 2.13 * u.degree




spectrum2 = Powerlaw()

shape2    = Disk_on_sphere()


source2   = ExtendedSource("gammaCygni",
                          spatial_shape=shape2,
                          spectral_shape=spectrum2)

fluxUnit = 1. / (u.TeV * u.cm**2 * u.s)


# Set spectral parameters (do it after the source definition to make sure
# the units are handled correctly)
spectrum2.K = 3.9e-12 * fluxUnit
spectrum2.K.fix = True

spectrum2.piv = 1.1 * u.TeV
spectrum2.piv.fix = True

spectrum2.index = -2.95
spectrum2.index.fix = True

shape2.lon0 = 305.27 * u.degree
shape2.lon0.fix = True

shape2.lat0 = 40.52 * u.degree
shape2.lat0.fix = True

shape2.radius = 0.63 * u.degree
shape2.radius.fix = True
shape2.radius.max_value = 0.63 * u.degree

lm = Model(source, source2, source3)
bin_list  =  "1c 1d 1e 1f 2c 2d 2e 2f 3c 3d 3e 3f 4c 4d 4e 4f 4g 5e 5f 5g 5h 6e 6f 6g 6h 7f 7g 7h 7i 8g 8h 8i 8j 9g 9h 9i 9j 9k 9l".split()


ra, dec = 307.17, 41.17
#data_radius = 6.0
model_radius = 8.0


fits_roi = HealpixMapROI(ra = ra, dec = dec, model_radius=model_radius, roifile="roi.fits")
hawc = HAL("HAWC", maptree, response,fits_roi)
hawc.set_active_measurements(bin_list=bin_list)
#hawc.display()

# Double check the free parameters
print("Likelihood model:\n")
print(lm)


# Set up the likelihood and run the fit
print("Performing likelihood fit...\n")
datalist = DataList(hawc)
jl = JointLikelihood(lm, datalist, verbose=True)
jl.set_minimizer("ROOT")
param_df, like_df = jl.fit()


TS = jl.compute_TS("veritasellipse", like_df)
print("TS:\n")
print(TS)

TS = jl.compute_TS("HAWCcocoon", like_df)
print("TS:\n")
print(TS) 


TS = jl.compute_TS("gammaCygni", like_df)
print("TS:\n")
print(TS)

contours = jl.get_contours(spectrum3.upper_bound_2,1e9,500e9,10)
plt.savefig('cocoon_up.pdf')
