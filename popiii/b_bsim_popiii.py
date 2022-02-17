import yt
from yt.units import dimensions
import matplotlib
matplotlib.use('ps')
from yt.mods import *
import numpy as na
import pylab as pl
#from string import rstrip
import fields_bfield
import tracer_def

def find_b(a, xi, bmag, phi_ff, mach, vt, xi_1u, chi):

	xi_1 = pow(pow(xi_1u,-0.5) + pow(xi,-0.5), -2.0)
	b1 = 1.7e-8 * pow(xi_1, 2.0/3.0)
	eps_b1 = 5.15e6 * pow(xi_1, 1.0/3.0)
	xi_rat = pow(xi_1 / xi, a)

	#print('xi_1 / xi_1u = ', xi_1/ xi_1u)

	A12 = 1.0 + (chi/(1.0/16.0)) * 8.12e-3 * (phi_ff * mach / a) * ( 1.0 - xi_rat) * (0.5 * vt * vt / eps_b1)
	B = b1 * np.sqrt(A12) * pow(xi/xi_1, (1.0+a)/2.0)

	return A12, B

def find_b_alt(a, xi, bmag, phi_ff, mach, vt, xi_1u, chi):
	chi_fac = .247
	if chi > (1.0/16.0):
		chi_fac = .312 
	xi_8 = xi / 1.e8
	g = pow(( 1.0 + 1.0 / pow(xi_8,0.5) ),2)
	#A12 = 1.0 + chi_fac * pow(g,1.0/3.0) * [ 1.0 - 1.0 / pow(g * xi_8, a) ] / a

	A12 = 1.0 + 0.25 / a * (chi/(1.0/16.0)) * pow(1.e8 / xi_1u, 1.0/3.0) 

	xi_1 = pow(pow(xi_1u,-0.5) + pow(xi,-0.5), -2.0)
	b1 = 1.7e-8 * pow(xi_1, 2.0/3.0)
	B = b1 * np.sqrt(A12) * pow(xi/xi_1, (1.0+a)/2.0)

	return A12, B


xi = [1.e8,2.96E+08, 4.18E+08, 5.89E+08, 8.30E+08, 1.17E+09, 1.65E+09, 2.32E+09, 3.28E+09, 4.62E+09,6.51E+09, 9.18E+09, 1.29E+10, 1.82E+10, 2.57E+10, 3.62E+10, 5.11E+10, 7.20E+10, 1.01E+11, 1.43E+11, 2.02E+11, 2.84E+11, 4.00E+11, 5.64E+11, 7.96E+11, 1.e25, 1.e30]

bmag = [1e-3, 7.20E-03, 8.59E-03, 1.05E-02, 1.29E-02, 1.61E-02, 2.03E-02, 2.53E-02, 3.05E-02, 3.74E-02, 4.66E-02, 5.83E-02, 7.39E-02, 9.20E-02, 1.13E-01, 1.42E-01, 1.78E-01, 2.16E-01, 2.59E-01, 3.14E-01, 3.79E-01, 4.52E-01, 5.54E-01, 6.99E-01, 8.09E-01, 8.09E-01,8.09E-01]

#xi = [1.e9, 1.e10, 1.e11]
#bmag = [ 8.59E-03, 1.05E-02, 1.29E-02]

vt = 2.e5

mach = 0.77 

xi = np.asarray(xi)
#bmag = np.asarray(bmag)
#vt = np.asarray(vt)
#mach = np.asarry(mach)

phi_ff = 4.7

xi_1u = 1.e8
chi = 1.0/16.0

a = 1.0/3.0
a12, bmag_theory_1 = find_b(a, xi, bmag, phi_ff, mach, vt, xi_1u, chi) 

a = 3.0/19.0
a12, bmag_theory_2 = find_b(a, xi, bmag, phi_ff, mach, vt, xi_1u, chi)

a = 4.0/57.0
a12_3, bmag_theory_3 = find_b(a, xi, bmag, phi_ff, mach, vt, xi_1u, chi)
a12_3_alt, bmag_alt_3 = find_b_alt(a, xi, bmag, phi_ff, mach, vt, xi_1u, chi)

a = 4.0/57.0
chi = 3.0/38.0
a12_4, bmag_theory_4 = find_b(a, xi, bmag, phi_ff, mach, vt, xi_1u, chi)
a12_4_alt, bmag_alt_4 = find_b_alt(a, xi, bmag, phi_ff, mach, vt, xi_1u, chi)

print('a12_3 = ', a12_3)
print('a12_3_alt = ', a12_3_alt)
print('a12_4 = ', a12_4)
print('a12_4_alt = ', a12_4_alt)

print('bmag_theory_3 =', bmag_theory_3)
print('bmag_alt_3 = ', bmag_alt_3)
print('bmag_theory_4 =', bmag_theory_4)
print('bmag_alt_4 = ', bmag_alt_4) 


print('min 1 = ', min(bmag_theory_1/bmag), 'max = ', max(bmag_theory_1/bmag))
print('min 2 = ', min(bmag_theory_2/bmag), 'max = ', max(bmag_theory_2/bmag))
print('min 3 = ', min(bmag_theory_3/bmag), 'max = ', max(bmag_theory_3/bmag))
print('min 4 = ', min(bmag_theory_4/bmag), 'max = ', max(bmag_theory_4/bmag))

fsize = 14

nmin = 1.e8
nmax = 1.e12

#bmin = 1e-5
#bmax = 1e1

bmin = 0.5
bmax = 2.5

#pl.plot(xi, bmag,'k:', linewidth=1.5, label='bmag')
pl.plot(xi, bmag_theory_1/bmag,'k', linewidth=1.5, label='a = 1/3 (flux-freezing)')
pl.plot(xi, bmag_theory_2/bmag,'b:', linewidth=3.0, label='a = 3/19 (homologous)')
pl.plot(xi, bmag_theory_3/bmag,'r--', linewidth=3.0, label='a = 4/57 (Jeans-regulated)')
pl.plot(xi, bmag_theory_4/bmag,'g-.', linewidth=3.0, label=r'a = 4/57 (Jeans-regulated, $\chi$=3/38)')
ax = pl.gca()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel(r'$\xi$', fontsize=fsize)
ax.set_ylabel(r'B / B$_{\rm sim}$', fontsize=fsize)
pl.xticks(fontsize=fsize)
pl.yticks(fontsize=fsize)
pl.axis((nmin, nmax, bmin, bmax))
pl.legend(loc=2, fontsize=10)
pl.savefig('b_bsim.eps', bbox_inches='tight')

