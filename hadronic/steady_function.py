import exceptions
import math

import astropy.units as u
import numpy as np
import warnings

from astromodels.core.units import get_units
from astromodels.functions.function import Function1D, FunctionMeta, ModelAssertionViolation
import warnings


from steady_description import *


#DiffCoeff_Rref = 1E9 # reference rigidity in Volts

class cutoff_pp(Function1D):
    r"""
    description                : transientModel paramters
    parameters                 :


        Q_inj_index :
                                desc: none 
                                initial value : 1.85
                                min : 0
                                max : 10
   
        E_inj_max :
                                desc: none 
                                initial value : 60e12
                                min : 1e10
                                max : 1e20



        k :
                                desc: none 
                                initial value : 1/30
                                is_normalization : True
                                min : None
                                max : None
                                delta: 0.001



        """

    __metaclass__ = FunctionMeta

    def _set_units(self, x_unit, y_unit):

            # This function can only be used as a spectrum,
            # so let's check that x_unit is a energy and y_unit is
            # differential flux

        if hasattr(x_unit, "physical_type") and x_unit.physical_type == 'energy':

            # Now check that y is a differential flux
            current_units = get_units()
            #should_be_unitless = y_unit * (current_units.energy * current_units.time * current_units.area)

        # For testing, I am putting all  dimensionless
        self.Q_inj_index.unit = u.dimensionless_unscaled
        self.E_inj_max.unit = u.dimensionless_unscaled
        self.k.unit = u.dimensionless_unscaled

    #def getModels(self, x, t_0, t_burst, dt_burst, DiffCoeff_0, DiffCoeff_index, Q_inj_index, E_inj_max):
        #protons=dN_dE_t0_p(x, t_0, DiffCoeff_0, DiffCoeff_index, Q_inj_index, E_inj_max)
        #return protons        
    # x in unit of KeV, change it to eV and convert flux in terms of keV-1cm-2s-1
    def evaluate( self, x, Q_inj_index, E_inj_max, k):
        x_eV = x * 1e3 
        gammas=phi_gamma_steady_p(x_eV, k, Q_inj_index, E_inj_max)
        gamma_flux = gammas/x_eV**2
        return gamma_flux * 1e3
