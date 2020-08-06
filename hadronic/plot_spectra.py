#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt

mpl.rc("font", family="serif", size=14)
    
from burst_description import*
from burst_function import*
from steady_description import*
from steady_function import*
    
def go():
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))


#ARGO flux points from plotdigitizer
    x = np.array ([0.437, 0.686, 0.991, 1.54, 2.86, 6.6])
    y = np.array([4.72e-11, 5.23e-11, 4.11e-11, 4.16e-11, 2.19e-11, 2.32e-11])* 0.624151
    y1 = np.array([1.68e-11, 3.72e-11, 2.76e-11, 2.88e-11, 9.16e-12, 1.19e-11])* 0.624151
    yerr1 = y-y1
    y2 = np.array([7.82e-11, 6.82e-11, 5.37e-11, 5.28e-11, 3.42e-11, 2.32e-11])* 0.624151
    yerr2 = y2 - y
    uplims = np.array([0, 0, 0, 0, 0, 2.32e-11])* 0.624151
    ax.errorbar(x, y, yerr=[yerr1,yerr2], uplims=uplims, fmt='s', color="k", label = "ARGO J2031+4157, 2014" )




#Data points from 2011 paper, Binita via slack.
    x = np.power(10, np.array( [9.14, 9.41, 9.65, 9.89, 10.26, 10.77]))/1e12
    y = np.power(10, np.array(  [-4.24832, -4.22174, -4.19497, -4.31060, -4.39871, -4.48776]) )*1e-12*1e6
    yerr = np.power( 10, np.array( [-4.13439, -4.10762, -4.13819, -4.19667, -4.31326, -4.34535])  )*1e-6 - y
    
    ax.errorbar(x, y, yerr=yerr, fmt='^', color="grey", label = "Fermi-LAT Cocoon, 2011" )


#4FGL
    LAT_4FGL_E_Left = np.array([50, 100, 300, 1e3, 3e3, 10e3, 30e3])*1e-6 # TeV

    LAT_4FGL_E_Right = np.array([100, 300, 1e3, 3e3, 10e3, 30e3, 300e3])* 1e-6 # TeV

    x = np.sqrt(LAT_4FGL_E_Left*LAT_4FGL_E_Right)
    xerr1 = x-LAT_4FGL_E_Left
    xerr2 = LAT_4FGL_E_Right-x
    LAT_4FGL_F = np.array([7.79840220e-11, 1.22992172e-10, 2.09725307e-10, 2.03366546e-10, 1.62071953e-10, 1.59847413e-10, 1.12895304e-10])*0.624151 # TeV cm^-2 s^-1

    LAT_4FGL_F_err = np.array([2.08776601e-11, 1.09858243e-11, 6.20437553e-12, 5.48733379e-12, 6.28723778e-12, 8.90127572e-12, 1.22363327e-11])* 0.624151 # TeV cm^-2 s^-1 

    ax.errorbar(x, y=LAT_4FGL_F, xerr=[xerr1, xerr2], yerr=[LAT_4FGL_F_err,LAT_4FGL_F_err], fmt='s', color="red", label = "Fermi-LAT Cocoon, 2018" )


#HAWC flux points
    x = np.array ([0.83, 1.20, 2.17, 4.49, 9.08, 17.72, 33.84, 60.90, 108.05, 191.7])
    y = np.array([4.1663330681e-11, 4.24957380986e-11, 2.45164387806e-11, 1.65417554382e-11, 1.11230854042e-11, 5.0911581727e-12, 2.50458962395e-12, 3.44247649348e-12, 3.27997251097e-12, 3.21645011935e-12])
    yerr1 = np.array([5.74938072e-12, 5.58512438e-12, 4.73626384e-12, 3.22141005e-12, 1.95680256e-12, 1.50178102e-12, 1.09073808e-12, 3.09822884e-12, 2.95197526e-12, 1.48494254e-12])
    yerr2 = np.array([5.84774935e-12, 5.59176155e-12, 4.86076856e-12, 3.10779800e-12, 2.02388102e-12, 1.46528007e-12, 1.14520887e-12, 0.00000000e+00, 0.00000000e+00, 1.53219434e-12])
    uplims = np.array([0, 0, 0, 0, 0, 0, 0, 3.09822884e-12, 2.95197526e-12, 0])
    #ax.errorbar(x, y, yerr=[yerr1,yerr2], fmt='o', uplims=uplims, color="green", label = "cocoon fluxpoints" )

#HAWC flux points last four bins combined
    x = np.array ([0.83, 1.20, 2.17, 4.49, 9.08, 17.72, 39.77, 122.12])
    y = np.array([4.18277121784e-11, 4.2446368961e-11, 2.46650166137e-11, 1.66425262843e-11, 1.12330927764e-11, 5.08165079672e-12, 2.41594514232e-12, 2.29373115472e-12])
    yerr1 = np.array([5.95078620e-12, 5.62583392e-12, 4.79169484e-12, 3.13839877e-12, 1.92521695e-12, 1.47463315e-12, 8.47518241e-13, 9.00180031e-13])
    yerr2 = np.array([5.76591186e-12, 5.84337101e-12, 4.64934865e-12, 3.09952335e-12, 1.90022128e-12, 1.45486529e-12, 8.43548285e-13, 9.10930504e-13])
    ax.errorbar(x, y, yerr=[yerr1,yerr2], fmt='o', color="blue", label = "HAWC Cocoon fluxpoints, this study" )

#HAWC flux points with flat bkg
    x = np.array ([0.82, 1.18, 2.14, 4.45, 9.03, 17.66, 33.76, 60.78, 107.84, 191.3])
    y = np.array([3.86099078887e-11, 3.97600687797e-11, 2.24692262751e-11, 1.49724273619e-11, 9.98894847244e-12, 4.25216484738e-12, 3.39071282953e-12, 3.05302718964e-12, 2.99918440875e-12, 3.03036276201e-12])
    yerr1 = np.array([5.97299958e-12, 5.76505119e-12, 4.93950285e-12, 3.18189470e-12, 1.88532841e-12, 1.46221495e-12, 3.05164155e-12, 2.74772447e-12, 2.69926597e-12, 1.49063540e-12])
    yerr2 = np.array([6.05160220e-12, 5.90747860e-12, 4.76272843e-12, 3.18513843e-12, 1.96750556e-12, 1.51401151e-12, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.44290013e-12])
    #ax.errorbar(x, y, yerr=[yerr1,yerr2], fmt='o', color="red", label = "cocoon fluxpoints" )


#HAWC flux points with flat_bkg last bins combined into half decades (ij and kl)
    x = np.array ([0.82, 1.18, 2.14, 4.45, 9.03, 17.66, 39.62, 121.75])
    y = np.array([3.86639696887e-11, 3.97477162228e-11, 2.24770582813e-11, 1.50091725874e-11, 1.00294402214e-11, 4.29453631271e-12, 1.87906462977e-12, 2.04907781681e-12])
    yerr1 = np.array([5.67651862e-12, 5.54774166e-12, 4.78472776e-12, 3.15360032e-12, 1.92040369e-12, 1.46505859e-12, 8.68875015e-13, 8.87655881e-13])
    yerr2 = np.array([5.80238824e-12, 5.75333872e-12, 4.87403982e-12, 3.11634030e-12, 1.87884764e-12, 1.44211542e-12, 8.86525138e-13, 8.99355721e-13])
    #ax.errorbar(x, y, yerr=[yerr1,yerr2], fmt='o', color="green", label = "HAWCcocoon fluxpoints" )



# Steady source scenario
    Q_inj_index = 2.127
    E_inj_max = 30e15
    k = 0.026

# gamma-ray energies
    x = np.logspace(7.5, 17, 400)*1e-3
# gamma-ray flux in steady model
    s = break_pp()
    y_Flux = s.evaluate(x, Q_inj_index, E_inj_max, k)
    y = y_Flux * x**2
    print x * 1e3
    print y* 1e3
    ax.errorbar(x/1e9, y/1e9, color="grey", label = "Burst Source Scenario" )


# Steady source scenario
    Q_inj_index_steady = 2.00
    E_inj_max_steady = 3e14
    k_steady = 0.00019

# gamma-ray energies
    x = np.logspace(7.5, 17, 400)*1e-3
# gamma-ray flux in steady model
    s = cutoff_pp()
    y_Flux_steady = s.evaluate(x, Q_inj_index_steady, E_inj_max_steady, k_steady)
    y = y_Flux_steady * x**2
    print x * 1e3
    print y* 1e3
    ax.errorbar(x/1e9, y/1e9, color="grey", fmt="--", label = "Steady source scenario" )

    ynew = s.evaluate(x, 1.81, 49e12, 0.00007)
    y1 = ynew * x**2
    print x * 1e3
    print y* 1e3
    #ax.errorbar(x/1e9, y1/1e9, color="red", label = "3ml fit Steady Source Scenario" )
#model
    #x1 = np.logspace(7.5, 17, 400)
    #y1 = phi_gamma_burst_p(x1, k, Q_inj_index, E_inj_max)
    #print y1
    #y = phi_gamma_burst_p(x_Energy, t_burst, dt_burst, DiffCoeff_0, DiffCoeff_index, Q_inj_index, E_inj_max)
    #ax.errorbar(x1/1e12, y1/1e12, color="black", label = "Ke's burst model" )  
# gamma-ray flux in steady model
    #y_Flux_steady = phi_gamma_steady_p(x_Energy, t_cocoon, DiffCoeff_0_steady, DiffCoeff_index_steady, Q_inj_index_steady, E_inj_max_steady)
    #ax.errorbar(x_Energy/1e12, y_Flux_steady/1e12, color = "yellow", label = "cutoff=112TeV" )
    #y1 = dN_dE_t0_p(x1, t0, Q_inj_index, E_inj_max)* x1**2
    #ax.errorbar(x1/1e12, y1/1e12, color="red", label = "proton burst" )  

    plt.yscale("log")
    plt.xscale("log")
    ax.set_xlim(0.00001, 316)
    ax.set_ylim( 1e-13, 1e-9)

    
    ax.set_xlabel("Energy [eV] ")
    ax.set_ylabel("E^2 dN/(dE dA dt) [TeV/(cm2 s)]")
    #ax.set_title("Preliminary Cocoon (V)HE gamma-ray spectra")
    L=ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig("both_spectra.png"  )
    fig.savefig("both_spectra.pdf"  )





if __name__=="__main__":

    go()
