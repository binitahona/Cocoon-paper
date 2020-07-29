import threeML
import astromodels
import numpy as np
import astropy.units as u
from astropy.constants import c

try:
    import naima
except ImportError:
    warnings.warn("The naima package is not available. Models that depend on it will not be available",
                  NaimaNotAvailable)
    has_naima = False
else:
    has_naima = True

from naima_functions import *
from threeML.io.progress_bar import progress_bar


def plot_spectrum( results, color, subplot, min, max, unit, sources=[], x=None, y=None, yerr=None, plot_style={} ):
    print sources
    print results.get_data_frame()
    threeML.plot_point_source_spectra( results, ene_min=min, ene_max=max, flux_unit=unit, grid=True, subplot=subplot, contour_colors=color, fit_colors=color, sources_to_use=sources, include_extended=True, num_ene=100, plot_style_kwargs=plot_style)
    if x is not None:
            subplot.errorbar(x, y, yerr=yerr, fmt='o', color=color, capsize=10)

def plot_model_spectrum( model_file, uncertainties, color, subplot, min, max, unit, sources=[], x=None, y=None, yerr=None, model=None, plot_style={} ):
    if model is None:
      model = astromodels.load_model( model_file )
    cov_list = []
    if not uncertainties:
      uncertainties=[]
      for (k, p) in model.free_parameters.items():
        uncertainties.append( p.delta*p.unit )
    
    for unc, (k, p)  in zip(uncertainties, model.free_parameters.items()):
      if p.unit != u.dimensionless_unscaled:
        unc = unc.to(p.unit).value
      if not p.has_transformation():
        cov_list.append(unc**2)
      else:  
        trafo=p.transformation.forward
        cov_list.append( (trafo( p.value + unc ) - trafo( p.value ))*(trafo( p.value )-trafo( p.value - unc )) )
    print cov_list    
    cov_matrix=np.diag( np.array(cov_list) )
    results = threeML.analysis_results.MLEResults(model, cov_matrix, {})  
    plot_spectrum( results, color, subplot, min, max, unit, sources, x, y, yerr, plot_style)


def plot_fitted_spectrum( fitsfile, color, subplot, min, max, unit, sources=[], x=None, y=None, yerr=None, plot_style={}  ):
    results = threeML.load_analysis_results( fitsfile )
    print sources
    plot_spectrum( results, color, subplot, min, max, unit, sources, x, y, yerr, plot_style)

def plot_fermipy_spectrum( npyfile, color, subplot, min, max, unit, sourcename, plotFluxPoints = True, plot_style={} , xunit = u.TeV ):
  c = np.load(npyfile).flat[0]
  
  if c["SpectrumType"] != "PowerLaw":
    print "plot_fermipy_results not implemented for SpectrumType", c["SpectrumType"]
  
  if plotFluxPoints:
    x=(c['e_ref']*u.MeV).to(xunit).value
    y=(c['e2dnde']*u.MeV/u.cm**2/u.s).to( unit ).value
    yerr=(c['e2dnde_err']*u.MeV/u.cm**2/u.s).to( unit ).value
    print x
    print y
    print yerr
  else:
    x=None
    y=None
    yerr=None
  
  spec = Powerlaw_linNorm()
  source = astromodels.PointSource( sourcename, ra=20., dec=20., spectral_shape = spec )
  spec.K = c['param_values'][0]/(u.MeV*u.cm**2*u.s)
  spec.K.bounds = ( 1e-50/(u.MeV*u.cm**2*u.s), 1/(u.MeV*u.cm**2*u.s))
  spec.index  = c['param_values'][1]
  spec.piv = c['param_values'][2]*u.MeV
  model = astromodels.Model(source)
  

  #trafo=spec.K.transformation.forward
  #unc_K=np.sqrt( ((c['param_covariance'][0,0] /(u.MeV**2*u.cm**4*u.s**2) ).to( u.keV**-2*u.cm**-4*u.s**-2 )).value )
  #K = spec.K.value
  #cov_KK =  (trafo( K + unc_K ) - trafo( K ))*(trafo( K )-trafo( K - unc_K )) 
  #cov_ii = c['param_covariance'][1,1]
  #cov_Ki = (np.sqrt( cov_KK ) / np.sqrt( unc_K ) *  c['param_covariance'][1,0] / (u.MeV*u.cm**2*u.s) ).to( u.keV**-1*u.cm**-2*u.s**-1 ).value
  
  #cov_matrix= np.array([ [ cov_KK, cov_Ki], [ cov_Ki, cov_ii] ] )
  
  cov_matrix = np.array( c["param_covariance"][0:2,0:2] )
  
  print cov_matrix
  
  cov_matrix[0,0]= (cov_matrix[0,0] / (u.MeV*u.cm**2*u.s)**2 ).to( (u.keV*u.cm**2*u.s)**-2 ).value
  cov_matrix[1,0]= (cov_matrix[1,0] / (u.MeV*u.cm**2*u.s) ).to( (u.keV*u.cm**2*u.s)**-1 ).value
  cov_matrix[0,1]= (cov_matrix[0,1] / (u.MeV*u.cm**2*u.s) ).to( (u.keV*u.cm**2*u.s)**-1 ).value
  
  print cov_matrix, spec.K

  
  results = threeML.analysis_results.MLEResults(model, cov_matrix, {})  
  plot_spectrum( results, color, subplot, min, max, unit, [sourcename], x, y, yerr, plot_style)

def plot_naima_ic_spectrum( color, subplot, min, max, unit, distance, We, alpha_1, alpha_2, Eb, Ec, B= 20*u.uG, n0=1/u.cm**3 ):
  energies=np.logspace(np.log10(min.to(u.keV).value), np.log10(max.to(u.keV).value), num=500)*u.keV

  BPL = naima.models.ExponentialCutoffBrokenPowerLaw(amplitude = 1/u.eV, e_0 = 1*u.GeV, e_break = Eb , alpha_1 =alpha_1, alpha_2 = alpha_2, e_cutoff = Ec, beta=1.0 )
  #PL = naima.models.PowerLaw(5e43*u.Unit('1/eV'), 1*u.GeV, 3.1)

  SYN = naima.models.Synchrotron(BPL, B=B)
  SYN.set_We(We)

  BREMS= naima.models.Bremsstrahlung(BPL, n0=n0)
  BREMS.set_We(We) 

  # Define energy array for synchrotron seed photon field and compute
  # Synchroton luminosity by setting distance to 0. The energy range should
  # capture most of the synchrotron output.
  #Lsy = SYN.flux(energies, distance=0*u.cm)
  #R = 15 * u.pc
  #phn_sy = Lsy / (4 * np.pi * R**2 * c) * 2.24

  #print c
  #print phn_sy[0]

  radiation=['CMB', 'FIR', 'NIR']#, ['UV', 30000*u.K, 1*u.eV*u.cm**-3] , ['SSC', energies, phn_sy]]

  IC = naima.models.InverseCompton(BPL, seed_photon_fields=radiation)
  IC.set_We( We) 
  
  fluxesS = SYN.sed(energies, distance=distance)
  fluxesIC = IC.sed(energies, distance=distance)
  fluxesB = BREMS.sed(energies, distance=distance)
  fluxesB[0:100]=np.zeros(100)

  subplot.plot( energies, (fluxesS+fluxesIC+fluxesB).to(unit).value, color, lw = 2, label= "Total emission")

  subplot.plot(energies, fluxesS.to(unit).value, color, lw=1, ls = "--", label="Synchroton" )
  
  subplot.plot(energies, fluxesB.to(unit).value, color, lw=1, ls = ":", label="Nonthermal Bremsstrahlung" )

  subplot.plot(energies, fluxesIC.to(unit).value, color, lw=1, label="IC (total)" )

  for ls, seed in enumerate(['CMB', 'FIR', 'NIR']):#, 'UV', 'SSC' ]):
    sed = (IC.sed(energies, seed=seed, distance=distance)).to(unit)
    subplot.plot(energies,sed,lw=1.5, ls="--", c=naima.plot.color_cycle[ls], label = seed)

    
def plot_naima_ebpl_spectrum( color, subplot, min, max, unit, distance, Wp, alpha_1, alpha_2, Eb, Ec, model="Pythia8" ):

    useLUT=True if model=="Pythia8" else False

    p = naima.models.ExponentialCutoffBrokenPowerLaw(amplitude = 1/u.eV, e_0 = 100*u.GeV, e_break = Eb , alpha_1 =alpha_1, alpha_2 = alpha_2, e_cutoff = Ec, beta=1.0 )
    PP = naima.models.PionDecay(p, nh=1.0 * u.cm** -3, nEpd=200, Epmin=1.0*u.GeV, Epmax=200*u.PeV, hiEmodel=model,useLUT=useLUT)
    PP.set_Wp( Wp )
    energies=np.logspace(np.log10(min.to(u.keV).value), np.log10(max.to(u.keV).value), num=500)*u.keV
    fluxes = energies*energies*PP.flux(energies, distance=distance)
    subplot.plot(energies, fluxes.to(unit).value, color, lw=2, label=model )

def plot_sampled_spectrum( fitsfile, color, subplot, min, max, unit, source, N=500, alpha=None ):
    results = threeML.load_analysis_results( fitsfile ) 
    best_model=results.optimized_model
    energies=np.logspace(np.log10(min.to(u.TeV).value), np.log10(max.to(u.TeV).value), num=50)*u.TeV

    the_source=best_model.sources[source].components['main']    
    fluxes_nu= energies*energies*the_source(energies.to(u.keV).value)/u.cm**2/u.s/u.keV
    
    
    n_samples = np.min( np.array([results.samples.shape[1], N]))

    samples = np.random.choice(range(results.samples.shape[1]), size=n_samples, replace=False, p=None)

    fluxes = np.zeros((n_samples, len(energies)) )

    with progress_bar(n_samples, title="Sampling") as p:

      for i in range(n_samples):
        for par, sample in zip( best_model.free_parameters, results.samples[:,samples[i]] ):
          best_model[par].value = sample
        the_source=best_model.sources[source].components['main']    
        fluxes[i,:] = (energies*energies*the_source(energies)).to(unit).value
        p.increase()

    if alpha is not None:
      subplot.plot(energies, fluxes.T, color=color, alpha=alpha, lw=1 )
    
    if alpha is None:
      flow = np.percentile( fluxes, 16.0, axis = 0 )
      fmed = np.percentile( fluxes, 50.0, axis = 0 )
      fhigh = np.percentile( fluxes, 84.0, axis = 0 )
    
      subplot.fill_between( energies, flow, fhigh, color=color, alpha = 0.5 )
    
      subplot.plot(energies, fmed, color=color, lw=1, ls=":" )
    
    subplot.plot(energies, fluxes_nu.to(unit).value, color=color, lw=2, label = source   )


def plot_random_spectra( fitsfile, color, subplot, min, max, unit, source, N=20, alpha = 0.01 ):
    results = threeML.load_analysis_results( fitsfile ) 
    best_model=results.optimized_model
    energies=np.logspace(np.log10(min.to(u.TeV).value), np.log10(max.to(u.TeV).value), num=50)*u.TeV

    the_source=best_model.sources[source].components['main']    
    fluxes_nu= energies*energies*the_source(energies.to(u.keV).value)/u.cm**2/u.s/u.keV
    
    
    n_samples = np.min( np.array([results.samples.shape[1], N]))

    samples = np.random.choice(range(results.samples.shape[1]), size=n_samples, replace=False, p=None)

    #fluxes = np.zeros((n_samples, len(energies)) )

    with progress_bar(n_samples, title="Sampling") as p:

      for i in range(n_samples):
        for par, sample in zip( best_model.free_parameters, results.samples[:,samples[i]] ):
          best_model[par].value = sample
        the_source=best_model.sources[source].components['main']    
        fluxes = (energies*energies*the_source(energies)).to(unit).value
        subplot.plot(energies, fluxes.to(unit).value, color=color, alpha=0.002 )
        p.increase()
    
    flow = np.percentile( fluxes, 16.0, axis = 0 )
    fmed = np.percentile( fluxes, 50.0, axis = 0 )
    fhigh = np.percentile( fluxes, 84.0, axis = 0 )
    
    
    subplot.fill_between( energies, flow, fhigh, color=color, alpha = 0.5 )
    
    subplot.plot(energies, fmed, color=color, lw=1, ls=":" )
    subplot.plot(energies, fluxes_nu.to(unit).value, color=color, lw=2, label = "Naima Fit!"   )


class Powerlaw_linNorm(astromodels.Function1D):
        r"""
        description :

            A simple power-law

        latex : $ K~\frac{x}{piv}^{index} $

        parameters :

            K :

                desc : Normalization (differential flux at the pivot value)
                initial value : 1.0
                is_normalization : True
                min : 1e-30
                max : 1e3
                delta : 0.1

            piv :

                desc : Pivot value
                initial value : 1
                fix : yes

            index :

                desc : Photon index
                initial value : -2
                min : -10
                max : 10

        tests :
            - { x : 10, function value: 0.01, tolerance: 1e-20}
            - { x : 100, function value: 0.0001, tolerance: 1e-20}

        """

        __metaclass__ = astromodels.FunctionMeta

        def _set_units(self, x_unit, y_unit):
            # The index is always dimensionless
            self.index.unit = u.dimensionless_unscaled

            # The pivot energy has always the same dimension as the x variable
            self.piv.unit = x_unit

            # The normalization has the same units as the y

            self.K.unit = y_unit

        # noinspection PyPep8Naming
        def evaluate(self, x, K, piv, index):

            xx = np.divide(x, piv)

            return K * np.power(xx, index)


class Log_parabola_linNorm(astromodels.Function1D):
    r"""
    description :

        A log-parabolic function. NOTE that we use the high-energy convention of using the natural log in place of the
        base-10 logarithm. This means that beta is a factor 1 / log10(e) larger than what returned by those software
        using the other convention.

    latex : $ K \left( \frac{x}{piv} \right)^{\alpha -\beta \log{\left( \frac{x}{piv} \right)}} $

    parameters :

        K :

            desc : Normalization
            initial value : 1.0
            is_normalization : True
            min : 1e-30
            max : 1e5

        piv :
            desc : Pivot (keep this fixed)
            initial value : 1
            fix : yes

        alpha :

            desc : index
            initial value : -2.0

        beta :

            desc : curvature (positive is concave, negative is convex)
            initial value : 1.0

    """

    __metaclass__ = astromodels.FunctionMeta

    def _set_units(self, x_unit, y_unit):

        # K has units of y

        self.K.unit = y_unit

        # piv has the same dimension as x
        self.piv.unit = x_unit

        # alpha and beta are dimensionless
        self.alpha.unit = u.dimensionless_unscaled
        self.beta.unit = u.dimensionless_unscaled

    def evaluate(self, x, K, piv, alpha, beta):

        #print("Receiving %s" % ([K, piv, alpha, beta]))

        xx = np.divide(x, piv)

        try:

            return K * xx ** (alpha - beta * np.log(xx))

        except ValueError:

            # The current version of astropy (1.1.x) has a bug for which quantities that have become
            # dimensionless because of a division (like xx here) are not recognized as such by the power
            # operator, which throws an exception: ValueError: Quantities and Units may only be raised to a scalar power
            # This is a quick fix, waiting for astropy 1.2 which will fix this

            xx = xx.to('')

            return K * xx ** (alpha - beta * np.log(xx))

    @property
    def peak_energy(self):
        """
        Returns the peak energy in the nuFnu spectrum

        :return: peak energy in keV
        """

        # Eq. 6 in Massaro et al. 2004
        # (http://adsabs.harvard.edu/abs/2004A%26A...413..489M)

        return self.piv.value * pow(10, ((2 + self.alpha.value) * np.log(10)) / (2 * self.beta.value))


class Cutoff_powerlaw_linE(astromodels.Function1D):
    r"""
    description :

        A power law multiplied by an exponential cutoff

    latex : $ K~\frac{x}{piv}^{index}~\exp{-x/xc} $

    parameters :

        K :

            desc : Normalization (differential flux at the pivot value)
            initial value : 1.0
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1

        piv :

            desc : Pivot value
            initial value : 1
            fix : yes

        index :

            desc : Photon index
            initial value : -2
            min : -10
            max : 10

        xc :

            desc : Cutoff energy
            initial value : 10.0

    """

    __metaclass__ = astromodels.FunctionMeta

    def _set_units(self, x_unit, y_unit):
        # The index is always dimensionless
        self.index.unit = u.dimensionless_unscaled

        # The pivot energy has always the same dimension as the x variable
        self.piv.unit = x_unit

        self.xc.unit = x_unit

        # The normalization has the same units as the y

        self.K.unit = y_unit

    # noinspection PyPep8Naming
    def evaluate(self, x, K, piv, index, xc):

        # Compute it in logarithm to avoid roundoff errors, then raise it
        log_v = index * np.log(x / piv) - x / xc

        return K * np.exp(log_v)

