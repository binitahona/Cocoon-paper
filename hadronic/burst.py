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


from burst_description import *
from burst_function import  *


# This disable momentarily the printing of warnings, which you might get
# if you don't have the Fermi ST or other software installed
with warnings.catch_warnings():

    warnings.simplefilter("ignore")

    from threeML import *

#Load data, detector response and roi
maptree = ... # maptree.hdf5 file
response = ... #  response.hdf5 file
roi = ..#Region of interest, roi.fits file

lm = load_model("hadronic/cocoon.yml")


fluxUnit = 1. / (u.TeV * u.cm**2 * u.s)
E = np.logspace(np.log10(1), np.log10(100), 100) * u.TeV

spectrumPP_BPL = break_pp()

shapeD = Gaussian_on_sphere()

hawcRA = 307.65
hawcDec = 40.93

from astromodels.core.sky_direction import SkyDirection
position = SkyDirection(hawcRA, hawcDec)

sourceCutP = ExtendedSource("cocoon_PP_BPL", spatial_shape=shapeD, spectral_shape=spectrumPP_BPL)


# NOTE: if you use units, you have to set up the values for the parameters

shapeD.lon0=hawcRA*u.degree
shapeD.lat0=hawcDec*u.degree
shapeD.lon0.fix=True
shapeD.lat0.fix=True

shapeD.sigma = 2*u.degree
shapeD.sigma.fix = True
shapeD.sigma.bounds = (0.2,2.0)*u.degree


s = spectrumPP_BPL


s.Q_inj_index = 2.14
s.Q_inj_index.fix = False
s.Q_inj_index.bounds = (0, 4)

s.E_inj_max = 1e18
s.E_inj_max.fix = True
#s.E_inj_max.bounds = (1e12, 1e23)

s.k = 0.03
s.k.fix = False
#s.k.bounds = (5, 200)

lm.add_source(sourceCutP)


bin_list  =  "1c 1d 1e 1f 2c 2d 2e 2f 3c 3d 3e 3f 4c 4d 4e 4f 4g 5e 5f 5g 5h 6e 6f 6g 6h 7f 7g 7h 7i 8g 8h 8i 8j 9g 9h 9i 9j 9k 9l".split()

ra, dec = 307.17, 41.17
#data_radius = 6.0
model_radius = 8.0


fits_roi = HealpixMapROI(ra = ra, dec = dec, model_radius=model_radius, roifile=roi)
hawc = HAL("HAWC", maptree, response,fits_roi)
hawc.set_active_measurements(bin_list=bin_list)
#hawc.display()


fluxkeV = u.cm**-2 / u.keV / u.s


x = [2.00e+06, 6.50e+06, 2.00e+07, 1.65e+08] *  u.keV
y = [(0.12693143, 0.10115737, 0.09976892, 0.07046372)] /(x*x) * fluxkeV
yerr = [(0.00342492, 0.00392419, 0.00555574, 0.00763732)] /(x*x)  * fluxkeV


#x = [2.00e+06, 6.50e+06, 2.00e+07, 1.65e+08] *  u.keV
#y = [(0.12693143, 0.10115737, 0.09976892, 0.07046372)] /(x*x) * fluxkeV
#yerr = [(0.00342492, 0.00392419, 0.00555574, 0.00763732)] /(x*x)  * fluxkeV

likeXY = UnresolvedExtendedXYLike("likeXY", x, y, yerr,  poisson_data=False, quiet=False, source_name="cocoon_PP_BPL")

# Double check the free parameters
print("Likelihood model:\n")
print(lm)


# Set up the likelihood and run the fit
print("Performing likelihood fit...\n")
datalist = DataList(hawc, likeXY)
jl = JointLikelihood(lm, datalist, verbose=True)
jl.set_minimizer("ROOT")
param_df, like_df = jl.fit()


jl.results.write_to("proton_bpl.fits", overwrite=True)
lm.save("proton_bpl.yml", overwrite=True)
