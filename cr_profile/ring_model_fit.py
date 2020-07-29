from hawc_hal import HAL, HealpixConeROI, HealpixMapROI
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc("font", family="serif", size=14)

import warnings

from threeML.plugin_prototype import PluginPrototype
from astromodels.core.sky_direction import SkyDirection
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from ring_function import *

# This disable momentarily the printing of warnings, which you might get
# if you don't have the Fermi ST or other software installed
with warnings.catch_warnings():

    warnings.simplefilter("ignore")

    from threeML import *


# Make sure that the HAWC plugin is available

assert is_plugin_available("HAWCLike"), "HAWCLike is not available. Check your configuration"

def go(args):

    new_model_reloaded = load_model("test_likelihood_0.6_free.yml") 

    #veritas spectrum
    spectrum = Cutoff_powerlaw()
    shape    = Gaussian_on_sphere()

    pw   = ExtendedSource("veritasellipse",
                              spatial_shape=shape,
                              spectral_shape=spectrum)

    fluxUnit = 1. / (u.TeV * u.cm**2 * u.s)

    # Set spectral parameters (do it after the source definition to make sure
    # the units are handled correctly)
    spectrum.K = 2.05e-13 * fluxUnit
    spectrum.K.fix = False
    spectrum.K.bounds = (1e-16 * fluxUnit, 1e-10*fluxUnit)

    spectrum.piv = 4.9 * u.TeV
    spectrum.piv.fix = True


    spectrum.index = -2.03
    spectrum.index.fix = False
    spectrum.index.bounds = (-4., 0.)

    spectrum.xc = 20 * u.TeV
    spectrum.xc.fix = False
    spectrum.xc.bounds = (1 * u.TeV, 400 * u.TeV)

    # Set spatial parameters (do it after the source definition to make sure
    # the units are handled correctly)
    shape.lon0 = 307.89 * u.degree
    shape.lon0.fix = True

    shape.lat0 = 41.58 * u.degree
    shape.lat0.fix = True

    shape.sigma = 0.27 * u.degree
    shape.sigma.fix = True
    shape.sigma.max_value = 0.27 * u.degree




    spectrum2 = Powerlaw()

    shape2    = Disk_on_sphere()


    gc   = ExtendedSource("gammaCygni",
                              spatial_shape=shape2,
                              spectral_shape=spectrum2)

    fluxUnit = 1. / (u.TeV * u.cm**2 * u.s)


    # Set spectral parameters (do it after the source definition to make sure
    # the units are handled correctly)
    spectrum2.K = 4.2e-12 * fluxUnit
    spectrum2.K.fix = True

    spectrum2.piv = 1.1 * u.TeV
    spectrum2.piv.fix = True

    spectrum2.index = -3.01
    spectrum2.index.fix = True

    shape2.lon0 = 305.27 * u.degree
    shape2.lon0.fix = True

    shape2.lat0 = 40.52 * u.degree
    shape2.lat0.fix = True

    shape2.radius = 0.63 * u.degree
    shape2.radius.fix = True
    shape2.radius.max_value = 0.63 * u.degree



    new_model_reloaded.add_source(gc)
    new_model_reloaded.add_source(pw)

    bin_list  =  "1c 1d 1e 1f 2c 2d 2e 2f 3c 3d 3e 3f 4c 4d 4e 4f 4g 5e 5f 5g 5h 6e 6f 6g 6h 7f 7g 7h 7i 8g 8h 8i 8j 9g 9h 9i 9j 9k 9l".split()

    #llh = HAWCLike("HAWC", args.mtfile, args.rsfile,fullsky=False)
    #llh.set_bin_list(bin_list)
    #llh.set_ROI(308.3,41.3,6,True)

    ra, dec = 307.17, 41.17
    #data_radius = 6.0
    model_radius = 8.0


    fits_roi = HealpixMapROI(ra = ra, dec = dec, model_radius=model_radius, roifile="/data/scratch/userspace/bhona/hal/new_roi.fits")
    hawc = HAL("HAWC", args.mtfile, args.rsfile,fits_roi)
    hawc.set_active_measurements(bin_list=bin_list)

  
 
    # Set up the likelihood and run the fit
    print("Performing likelihood fit...\n")
    datalist = DataList(hawc)
    jl = JointLikelihood(new_model_reloaded, datalist, verbose=True)
    jl.set_minimizer("ROOT")
    param_df, like_df = jl.fit()

    TS = jl.compute_TS("Ring0", like_df)
    print("TS:\n")
    print(TS)

    TS = jl.compute_TS("Ring1", like_df)
    print("TS:\n")
    print(TS)


    TS = jl.compute_TS("Ring2", like_df)
    print("TS:\n")
    print(TS)


    TS = jl.compute_TS("Ring3", like_df)
    print("TS:\n")
    print(TS)


    #jl.get_errors()
    jl.results.write_to("results/ringfit.fits", overwrite=True)
    #new_model_reloaded.save("results/ringfit.yml", overwrite=True)



if __name__=="__main__":

    import argparse

    p = argparse.ArgumentParser(description="Example spectral fit with LiFF")
    p.add_argument("-m", "--maptreefile", dest="mtfile",
                   help="LiFF MapTree ROOT file", required=True)
    p.add_argument("-r", "--responsefile", dest="rsfile",
                   help="LiFF detector response ROOT file", required=True)
    p.add_argument("-n", "--ntransit", dest="ntransit", default=1, type=float,
                   help="Number of transits to use")
    p.add_argument("--startBin", dest="startBin", default=1, type=int,
                   help="Starting analysis bin [0..9]")
    p.add_argument("--stopBin", dest="stopBin", default=9, type=int,
                   help="Stopping analysis bin [0..9]")
    p.add_argument("--RA", dest="RA", default=307.93, type=float,
                   help="Source RA in degrees (Hegra default)")
    p.add_argument("--Dec", dest="Dec", default=41.51, type=float,
                   help="Source Dec in degrees (Hegra default)")
    p.add_argument("--nMCMC", dest="nMCMC", default=0, type=int,
                   help="Number of MCMC samples to obtain (default=0).")
    p.add_argument("--ROI", dest="ROI", type=int)
    args = p.parse_args()

    go(args)

