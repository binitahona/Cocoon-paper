import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc("font", family="serif", size=14)

import warnings

from threeML.plugin_prototype import PluginPrototype
from astromodels.core.sky_direction import SkyDirection
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

#from ring_function import *

# This disable momentarily the printing of warnings, which you might get
# if you don't have the Fermi ST or other software installed
with warnings.catch_warnings():

    warnings.simplefilter("ignore")

    from threeML import *


# Make sure that the HAWC plugin is available

assert is_plugin_available("HAWCLike"), "HAWCLike is not available. Check your configuration"

def go(args):

    new_model_reloaded = load_model(args.model) 


    bin_list  =  "1c 1d 1e 1f 2c 2d 2e 2f 3c 3d 3e 3f 4c 4d 4e 4f 4g 5e 5f 5g 5h 6e 6f 6g 6h 7f 7g 7h 7i 8g 8h 8i 8j 9g 9h 9i 9j 9k 9l".split()

    llh = HAWCLike("HAWC", args.mtfile, args.rsfile,fullsky=True)
    llh.set_bin_list(bin_list)

    #llh.set_template_ROI("roi.fits", 0.5, True)

    llh.set_ROI(307.17, 41.17, 6, True)
    # Double check the free parameters
    print("Likelihood model:\n")
    print(new_model_reloaded)



    # Set up the likelihood and run the fit
    print("Performing likelihood fit...\n")
    datalist = DataList(llh)
    jl = JointLikelihood(new_model_reloaded, datalist, verbose=False)
    llh.calc_TS()



    llh.write_model_map(args.modelmapname,False)
    #llh.write_residual_map(args.residualmapname)





if __name__=="__main__":

    import argparse

    p = argparse.ArgumentParser(description="Example spectral fit with LiFF")
    p.add_argument("-m", "--maptreefile", dest="mtfile",
                   help="LiFF MapTree ROOT file", required=True)
    p.add_argument("-r", "--responsefile", dest="rsfile",
                   help="LiFF detector response ROOT file", required=True)
    p.add_argument("-yml", "--yml", dest="model",
                   help="yml model", required=True)
    p.add_argument("-mname", "--modelmapname", dest="modelmapname",
                   help="name of the model map root file", required=True)
    p.add_argument("-rname", "--residualmapname", dest="residualmapname",
                   help="name for residual map root file", required=True)

    args = p.parse_args()

    go(args)

