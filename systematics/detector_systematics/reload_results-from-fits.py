from threeML import *

import glob

#fits_file_name = ..  #name of the fits files from the results directory
for fits_file_name in glob.glob("results/*fits"):
    
    results_reloaded = load_analysis_results(fits_file_name)

    print fits_file_name
    results_reloaded.display()
    print "\n\n"
