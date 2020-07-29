import numpy as np
from scipy.special import erf
from scipy import integrate
from scipy.special import hyp1f1

# constants in cgs units
Myr = 3.1536e13
parsec = 3.0857e18
c_cgs = 2.9979E10
sigma_pp_inelastic_PeV = 59.2 * 1e-27 # at PeV based on eq 79 of Kelner (2006)
sigma_T = 6.6525E-25

# parameters in cgs units
R_cocoon = 60 * parsec
D_cocoon = 1.4e3 * parsec
n_gas = 30.
B_mag = 20e-6
Q_inj_0 = (3e39 / 1.6e-12) # to be normalized

b_cooling_b_p = n_gas * sigma_pp_inelastic_PeV * c_cgs
b_cooling_b_e = 4./3 * c_cgs * sigma_T * B_mag**2 / (8 * np.pi) / (0.511e6)**2 / 1.6e-12


# transientModel paramters
DiffCoeff_Rref = 1E9 # reference rigidity in Volts
DiffCoeff_index = 0.55
DiffCoeff_0 = 1e25 * DiffCoeff_Rref**(-DiffCoeff_index) #2e24
#tau = 0.8 * Myr
tau = 3 * Myr

# dE/dt = b E
def b_E_p(E):
    return  b_cooling_b_p * E

# diffusion distance for proton with E today traveling over t0
def lambda_distance_p(E, t0, DiffCoeff_0, DiffCoeff_index):
    # with b(E) = bE for hadronic interaction
    lambda_square = DiffCoeff_0 / b_cooling_b_p / DiffCoeff_index * E ** DiffCoeff_index * (np.exp(b_cooling_b_p * t0 * DiffCoeff_index) - 1)
    return lambda_square ** 0.5


# proton energy at t0
def E_0_p(E, t0):
    return E * np.exp(b_cooling_b_p * t0)



# proton spectrum in transient scenario
# dN_p/dE_p integrated over the cocoon volumn for burst at time t
# DiffCoeff = DdiffCoeff_0 E ^ DiffCoeff_index = 3x10^28 (E / 3GV)^ DiffCoeff_index for ISM
# Q_inj_index is proton injection index
# E_inj_max is proton maximum energy
def dN_dE_t0_p(E, t0, k, Q_inj_index, E_inj_max):
    lambda_0 = lambda_distance_p(E, t0, DiffCoeff_0, DiffCoeff_index)
    x = R_cocoon / lambda_0 / 2
    Fcc = erf(x) - 2 * x / np.pi**0.5 * np.exp(-x**2)
    E_0 = E * np.exp(b_cooling_b_p * t0)
    Q_inj = Q_inj_0 * k
    return E_0 / E * Q_inj * E_0**(-Q_inj_index) * np.exp(-E_0 / E_inj_max) * Fcc


def dN_dE_steady(E, dN_dE_t0, k, Q_inj_index, E_inj_max):
    logt, dlogt = np.linspace(-6, 0, 2049, retstep = True)
    t = 10.**logt * tau
    inte = E * 0.
    for i, _E in enumerate( E ):
        inte[i] = integrate.romb(y=dN_dE_t0(_E, t, k, Q_inj_index, E_inj_max) * t, dx=dlogt) * np.log(10.)
    return inte

# Gamma-ray flux from pp
# eq 79 of Kelner (2006)
def sigma_pp_inelastic(E):
    L = np.log(E / 1e12)
    Eth = 1.22e9
    thermalFactor = (1 - (Eth / E)**4)**2 * np.heaviside(E - Eth, 0.)
    return (34.3 + 1.88 * L + 0.25 * L**2) * 1e-27 * thermalFactor # cm^2


# eq 58 of 0606058 Kelner et al (2006)
def F_gamma(x, E_p):
    L = np.log(E_p / 1e12)
    
    B_g = 1.30 + 0.14 * L + 0.011 * L**2
    beta_g = 1. / (1.79 + 0.11 * L + 0.008 * L**2)
    k_g = 1. / (0.801 + 0.049 * L + 0.014 * L**2)
    
    y = x ** beta_g
    prefactor = B_g * np.log(x) / x * ((1 - y) / (1 + k_g * y * (1- y)))**4
    first_term = 1 / np.log(x)

    second_term = -4 * beta_g * y / (1 - y)
    third_term = - 4 * k_g * beta_g * y * (1 - 2 * y) / (1 + k_g * y * (1 - y))
    return prefactor * (first_term + second_term + third_term)



# gamma-ray spectrum in the steady scenario
def phi_gamma_steady_p(E_gamma, k, Q_inj_index, E_inj_max):
    phi_gamma_ = E_gamma * 0
    for i, E_gamma_ in enumerate(E_gamma):
        E_p_max = E_inj_max * 10. if E_gamma_ < E_inj_max else E_gamma_ * 10.
        logE_p, dlogE_p = np.linspace(np.log10(E_gamma_*(1+1e-10)), np.log10(E_p_max), 2049, retstep=True)
        E_p = 10.**logE_p
        phi_gamma_[i] = integrate.romb(y=dN_dE_steady(E_p, dN_dE_t0_p, k, Q_inj_index, E_inj_max) * sigma_pp_inelastic(E_p) * F_gamma(E_gamma_ / E_p, E_p), dx=dlogE_p) * c_cgs * n_gas * np.log(10.)
    return phi_gamma_ / (4 * np.pi * D_cocoon**2) * E_gamma**2
