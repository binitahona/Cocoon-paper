import numpy as np
from scipy.special import erf
from scipy import integrate
from scipy.special import hyp1f1

# constants in cgs units
Myr = 3.1536e13
parsec = 3.0857e18
c_cgs = 2.9979E10
sigma_T = 6.6525E-25

# parameters in cgs units
R_cocoon = 60 * parsec
D_cocoon = 1.4e3 * parsec
n_gas = 30.
B_mag = 20e-6
Q_inj_0 = (3e39 / 1.6e-12) # to be normalized

b_cooling_b_e = 4./3 * c_cgs * sigma_T * B_mag**2 / (8 * np.pi) / (0.511e6)**2 / 1.6e-12


# transientModel paramters
DiffCoeff_Rref = 1E9 # reference rigidity in Volts
DiffCoeff_index = 0.33
DiffCoeff_0 = 1e28 * DiffCoeff_Rref**(-DiffCoeff_index) #2e24
#tau = 0.8 * Myr
tau = 3 * Myr


def b_E_e(E):
    return b_cooling_b_e * E**2

def lambda_distance_e(E, t0, DiffCoeff_0, DiffCoeff_index):
    lambda_square = DiffCoeff_0 / b_cooling_b_e / (1 - DiffCoeff_index) * E ** (DiffCoeff_index - 1) * (1 - (1 - E * b_cooling_b_e * t0)**(1 - DiffCoeff_index))
    return lambda_square ** 0.5

def E_0_e(E, t0):
    return E / (1 - E * b_cooling_b_e * t0)

# electron; assuming E is one value, not array
def dN_dE_t0_e(E, t0, Q_inj_index, E_inj_max):
    dN_dE = t0 * 0.
    if isinstance(t0, (list, tuple, np.ndarray)):
        aliveEIndex = np.where(t0 < 1. / E / b_cooling_b_e)
        #if len(aliveEIndex) == 1:
            #print E, t0 / Myr
        lambda_0 = lambda_distance_e(E, t0[aliveEIndex], DiffCoeff_0, DiffCoeff_index)
        x = R_cocoon / lambda_0 / 2
        Fcc = erf(x) - 2 * x / np.pi**0.5 * np.exp(-x**2)
        E_0 = E_0_e(E, t0[aliveEIndex])
        dN_dE[aliveEIndex] = E**(-2) * Q_inj_0 * E_0 ** (2-Q_inj_index) * np.exp(-E_0 / E_inj_max) * Fcc

    else:
        if t0 < 1. / E / b_cooling_b_e:
            lambda_0 = lambda_distance_e(E, t0, DiffCoeff_0, DiffCoeff_index)
            x = R_cocoon / lambda_0 / 2
            Fcc = erf(x) - 2 * x / np.pi**0.5 * np.exp(-x**2)
            E_0 = E_0_e(E, t0)
            dN_dE = E**(-2) * Q_inj_0 * E_0 ** (2-Q_inj_index) * np.exp(-E_0 / E_inj_max) * Fcc
    return dN_dE


def dN_dE_steady(E, Q_inj_index, E_inj_max):
    logt, dlogt = np.linspace(-6, 0, 2049, retstep = True)
    t = 10.**logt * tau
    inte = E * 0.
    for i, _E in enumerate( E ):
        inte[i] = integrate.romb(y=dN_dE_t0_e(_E, t, Q_inj_index, E_inj_max) * t, dx=dlogt) * np.log(10.)
    return inte

