
from threeML import *

from ring_function import *

#fitsfile = raw_input()

results_reloaded = load_analysis_results("results/ringfit.fits")

results_reloaded.display()
for i in range(0,4):
	fluxes = results_reloaded.get_flux(1 * u.TeV, 200 * u.TeV, sources=("Ring%d" %i),include_extended=True)

	print fluxes

