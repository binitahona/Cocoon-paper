import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import ascii
import numpy as np
import math
import matplotlib as mpl
mpl.use("TkAgg")
from scipy.optimize import curve_fit

#unit for fluxes: erg / (cm2 s)
flux = np.array([2.3, 14, 14.8, 11])*1e-12

flux_l = np.array([1.7, 2.2, 2.7, 4.0])*1e-12

flux_h = np.array([7, 2.9, 4, 7])*1e-12

#energy range 1 tev to 200 tev
#flux = np.array([2.4, 14.2, 14.2, 11])*1e-12

#flux_l = np.array([1.8, 2.8, 3.2, 5.0])*1e-12

#flux_h = np.array([7, 3.3, 4, 8])*1e-12

flux_low = flux - flux_l

flux_high = flux + flux_h

#L=luminosity, given by d^2 4pi f, where d = 1.4 kpc for cygnus
d = 1.4 * 3.086 * 1e21

L = np. array(flux * d * d * 4 * math.pi) 

L_low = np. array(flux_low * d * d * 4 * math.pi)

L_high = np. array(flux_high * d * d * 4 * math.pi)
print L
print L - L_low
print L_high - L


#Unit of CR density is eV/cm3, formula from Aharonian SFR paper, eq. 7
L_f = np.array(L)/(10**34)

L_flow = np.array(L_low)/(10**34)

L_fhigh = np.array(L_high)/(10**34)

M_f = np.array([0.08, 0.24, 0.4, 0.33])**-1

n = ((1.5)/1.5)**-1

k = 1.8 * 10**-2 * n

w_cr = np.array(k * L_f * M_f)

werr1 = np.array((k * L_f * M_f)-(k * L_flow * M_f))

werr2 = np.array((k * L_fhigh * M_f)- (k * L_f * M_f))

print w_cr
print werr1
print werr2

r = np.array([0.6, 1.2, 1.8, 2.4])

r_0 = np.array([0, 0.6, 1.2, 1.8])
r_l = np.array(r*1.4*1000* math.pi)/180
r_0l = np.array(r_0*1.4*1000* math.pi)/180

print r_l
print r_0l

r_pc = (r_l + r_0l)/2

r_geo = np.array([7.33, 20.73, 35.91, 50.79])
xerr2 = r_l - r_geo
xerr1 = r_geo - r_0l

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.errorbar(r_geo, w_cr, yerr=[werr1,werr2], xerr=[xerr1,xerr2], fmt='o', color='g', label ='CR density (>10 TeV)')


gev_cr = np.array([0.13, 0.065, 0.041, 0.035])
gev_cr_err = np.array([0.013, 0.007, 0.006,0.004])
ax.errorbar(r_geo, gev_cr, yerr=[gev_cr_err,gev_cr_err], xerr=[xerr1,xerr2], fmt='^', color='k', label ='CR density (>100 GeV), Aharonian et al, 2019')

def const_func(x, n):
    return  (x**0) * n

def main():
    #a = 0.031
    #b = 10
    x = r_geo
    y = np.array([0.012138586933140495, 0.024629016965792314, 0.015621833618302551, 0.01407372398045275])
    y_h = np.array([0.008971999037538626, 0.003870274094624506, 0.0028499291060416817, 0.005117717811073728])
    y_l = np.array([0.03694352544868845, 0.00510172494291412, 0.004222117194135821, 0.008956006169379017])
    udata = (y_h + y_l)/2
    popt, pcov = curve_fit(const_func, x, y, sigma = udata, absolute_sigma=True)
    x_fit= np.arange(0.1, 60)
    y_fit = const_func(x_fit, *popt)
    plt.plot(x_fit, y_fit, label='Constant profile' )
    print (popt)
    print 'Best fit constant =', popt[0]
    print 'Uncertainty on constant =', np.sqrt(pcov[0, 0])

    chi = (y - const_func(x, *popt)) / udata
    chi2 = (chi ** 2).sum()
    dof = len(x) - len(popt)
    factor = (chi2 / dof)
    pcov_sigma = pcov / factor
    print chi2
    print dof
    print factor
    print np.sqrt(pcov_sigma[0, 0])

main()

def test_func(x,n):
    return (1/np.sqrt(60**2 - x**2) * np.log((60 + np.sqrt(60**2 - x**2))/x)) * n

def main():
    #a = 0.031
    #b = 10
    x = r_geo
    y = np.array([0.012138586933140495, 0.024629016965792314, 0.015621833618302551, 0.01407372398045275])
    y_h = np.array([0.008971999037538626, 0.003870274094624506, 0.0028499291060416817, 0.005117717811073728])
    y_l = np.array([0.03694352544868845, 0.00510172494291412, 0.004222117194135821, 0.008956006169379017])
    udata = (y_h + y_l)/2
    popt, pcov = curve_fit(test_func, x, y, sigma = udata, absolute_sigma=True)
    x_fit= np.arange(0.1, 60)
    y_fit = test_func(x_fit, *popt)
    plt.plot(x_fit, y_fit, label='1/r profile' )
    print (popt)
    print 'Best fit constant =', popt[0]
    print 'Uncertainty on constant =', np.sqrt(pcov[0, 0])

    chi = (y - const_func(x, *popt)) / udata
    chi2 = (chi ** 2).sum()
    dof = len(x) - len(popt)
    factor = (chi2 / dof)
    pcov_sigma = pcov / factor
    print chi2
    print dof
    print factor
    print np.sqrt(pcov_sigma[0, 0])
main()

# local cosmic ray density
x_constant = np.array([0, 7.33, 10, 12, 16, 20.73, 24, 28, 32, 35.91, 40, 45, 50.79, 60])
#y_constant = np.array([0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008])
y_constant = np.array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])
ax.errorbar(x_constant, y_constant, fmt='--', color='k', label ='local CR (>10 TeV)')

ax.set_ylabel("CR Energy density (eV/cm^3)")

ax.set_xlabel("Projected radius (pc)")
plt.rcParams.update({'font.size': 11})

plt.rcParams.update({'axes.titlesize': 16})
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

ax.set(xlim=(0, 60), ylim=(0, 0.15))
plt.legend()
plt.show()

fig.savefig("paper_cr_plot.png")
