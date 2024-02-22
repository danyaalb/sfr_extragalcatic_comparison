from __future__ import division
import pandas as pd
import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy.table import Table
import matplotlib.pyplot as plt

"""
This first block of code imports a table of the kinematic propeties of the CMZ catalogue.
Each structure in the catalogue has a number of kinematic properties that we have measured.
For our purposes here we only need a few of these, and so we use a Python package
called 'Pandas' to drop certain columns and rows in the table. Whenever we drop
rows or columns, it's best to re-index the table, so that the rows are indexed as
[0, 1, 2, 3, ...,]
"""

cmz_clouds_kinematic = Table.read('cmz_cloud_kinematic_properties.tex').to_pandas()
cmz_clouds_kinematic = cmz_clouds_kinematic[['Structure', 'HNCO mom1', 'HNCO fwhm']]
cmz_clouds_kinematic = cmz_clouds_kinematic.drop(labels=0, axis=0)
cmz_clouds_kinematic = cmz_clouds_kinematic.reset_index(drop=True)
# Some structures do not have valid kinematic measurements.
# We make a list of these structures, and drop them from our catalogue.
droplist = cmz_clouds_kinematic.index[cmz_clouds_kinematic['HNCO mom1'] == '-'].tolist()
cmz_clouds_kinematic = cmz_clouds_kinematic.drop(labels=droplist, axis=0)
cmz_clouds_kinematic = cmz_clouds_kinematic.reset_index(drop=True)


"""
This second block is similar to the above, but for the measured physical properties
such as mass, radius, etc. Here we have to drop the same entries as before, so that
the two tables are consistent with each other.
"""

cmz_clouds_physical = Table.read('cmz_cloud_physical_properties.tex').to_pandas()
cmz_clouds_physical = cmz_clouds_physical[['Structure', 'Area', 'l', 'b',
                        'Median column density', 'Mass', 'Radius',
                        'Median dust temperature']]
cmz_clouds_physical = cmz_clouds_physical.drop(labels=0, axis=0)
cmz_clouds_physical = cmz_clouds_physical.reset_index(drop=True)
cmz_clouds_physical = cmz_clouds_physical.drop(labels=droplist, axis=0)
cmz_clouds_physical = cmz_clouds_physical.reset_index(drop=True)

"""
Now we select specific columns from our tables to obtain our desired values of
Mass, Radius, etc. We use the Astropy units functionality to do the necessary unit
conversions.
"""

Mass     = cmz_clouds_physical['Mass'].values.astype(float)*u.Msun
Radius   = cmz_clouds_physical['Radius'].values.astype(float)*u.pc
Area     = cmz_clouds_physical['Area'].values.astype(float)*(u.pc*u.pc)
# sigma = velocity dispersion. Our measured value is actually the line width
# (full width at half maximum / FWHM), so we need to convert it as shown below.
sigma    = np.around((cmz_clouds_kinematic['HNCO fwhm'].values.astype(float)*(u.km/u.s))/ (2*np.sqrt(2*np.log(2))), 1)
surf_den = Mass.to(u.g) / Area.to(u.cm*u.cm)
sigma2_r = np.around((sigma**2 / Radius), 1)

# Plot surface density vs. sigma^2/r
plt.plot(np.log10(surf_den.value),np.log10(sigma2_r.value),ls='', marker='o', color='g', label='CMZ clouds', markersize='3', markeredgecolor='k')

"""
The below data are from the Galactic Ring Survey (GRS). The catalogue here contains
the equivalent measurements (surface density & sigma^2/r) towards clouds in the
disc of the Milky Way.
"""
grs = np.loadtxt('GRS.txt')
for i in range(0,len(grs),1):
     x = grs[i][0]
     y = grs[i][1]
     if i==0:
          plt.plot(x,y,ls='',marker=',',color='k',label='Heyer+ 09',markersize='4')
     else:
          plt.plot(x,y,ls='',marker=',',color='k',markersize='3')

"""
The below data are from the a paper by Faesi et al. (2018). The catalogue here
contains  the equivalent measurements (surface density & sigma^2/r) towards clouds
in an external galaxy, NGC 300.
"""
f18 = np.loadtxt('f18.txt')
for i in range(0,len(f18),1):
     x = f18[i][0]
     y = f18[i][1]
     if i==0:
          plt.plot(np.log10(((x*(u.Msun/(u.pc*u.pc))).to(u.g/(u.cm*u.cm))).value),np.log10(y),ls='',marker='s',color='b',label='Faesi+ 18',markersize='3')
     else:
          plt.plot(np.log10(((x*(u.Msun/(u.pc*u.pc))).to(u.g/(u.cm*u.cm))).value),np.log10(y),ls='',marker='s',color='b',markersize='3')

"""
I've commented this block out for now.
This will add the lines/curves for Virial equilibrium and pressure-bounded equilibrium.
We can leave this out for now and discuss later, but the feel free to un-comment
the below code and plot it, just to see what it looks like.
"""

"""
Gamma   = 0.73
n = 0
Y = np.zeros(shape=[5000])
Sig = np.zeros(shape=[5000])
Pe = 0

for r in np.logspace(-2.5,1.5,num=5000):
    Y[n] = (1/3)*((np.pi*Gamma*G*r*10) + ((4*(Pe*(k*1e6)))/(r*10)))
    Sig[n] = r
    n = n+1

plt.plot(np.log10(Sig), np.log10((Y*pc)/1e6), color='k', ls='--', alpha=0.4)

n = 0
Pe = [1e5, 1e6, 1e7, 1e8, 1e9]

for i in range(0,len(Pe),1):
    for r in np.logspace(-2.5,1.5,num=5000):
        Y[n] = (1/3)*((np.pi*Gamma*G*r*10) + ((4*(Pe[i]*(k*1e6)))/(r*10)))
        Sig[n] = r
        n = n+1

    plt.plot(np.log10(Sig), np.log10((Y*pc)/1e6), color='k', alpha=0.4)
    n=0
plt.text(-1.3, 2.2, '10$^{8}$ K cm$^{-3}$', color='k', fontsize=12)
plt.text(-1.7, 1.65, '10$^{7}$ K cm$^{-3}$', color='k', fontsize=12)
plt.text(-2.1, 1.0, '10$^{6}$ K cm$^{-3}$', color='k', fontsize=12)
"""

# Update some plotting parameters and save the figure.
plt.ylabel('log $\sigma$$^{2}$/ R (km$^{2}$ s$^{-2}$ pc$^{-1}$)')
plt.xlabel('log $\Sigma$ (g cm$^{-2}$)')
plt.legend(loc='lower right', numpoints=1, ncol=1, fontsize=10)
plt.ylim(-1.0, 3.0)
plt.xlim(-2.5, 0.5)
plt.savefig('CMZ_comparison_plot_test.pdf')
plt.clf()
