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

new_model_reloaded = load_model("models/pwn+cocoon+cygni.yml")
subtract = load_model("models/pwn+cygni.yml")

bin_list  =  "1c 1d 1e 1f 2c 2d 2e 2f 3c 3d 3e 3f 4c 4d 4e 4f 4g 5e 5f 5g 5h 6e 6f 6g 6h 7f 7g 7h 7i 8g 8h 8i 8j 9g 9h 9i 9j 9k 9l".split()


ra, dec = 307.17, 41.17
#data_radius = 6.0
model_radius = 8.0


fits_roi = HealpixMapROI(ra = ra, dec = dec, model_radius=model_radius, roifile="roi.fits")
hawc = HAL("HAWC", maptree, response,fits_roi)
hawc.set_active_measurements(bin_list=bin_list)
#hawc.display()

# Set up the likelihood and run the fit
print("Performing likelihood fit...\n")
datalist = DataList(hawc)
jl = JointLikelihood(new_model_reloaded, datalist, verbose=True)
jl.set_minimizer("ROOT")
param_df, like_df = jl.fit()

TS = jl.compute_TS("HAWCcocoon", like_df)
print("TS:\n")
print(TS) 


raplot = llh.plot_radial_profile(307.65, 40.93, max_radius=2.3)
raplot.savefig("pwn+cygni_subtracted.png")
