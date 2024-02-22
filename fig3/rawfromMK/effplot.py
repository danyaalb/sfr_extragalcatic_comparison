# Import libraries
from numpy import *
from matplotlib.pyplot import *
from astropy.io import ascii
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.colorbar as cbar
from interpshift import interpshift

# Define unit conversions
pc=3.09e18
kpc=1e3*pc
msun=1.99e33
myr=3.16e13
g=6.67e-8
msunpc2=msun/pc**2
msunpc2myr=msunpc2/myr

# Choice of CO alpha factor scalings NORMALIZED TO DADDI'S CONVENTIONS
alphasb = 1.0      # starburst value
alphahiz = 1.0     # high-z disk
alphaz0 = 1.0      # z=0 scaling
alphathings = 4.6/4.4 # THINGS data value

# Parameters
sigmagmc=85*msun/pc**2    # GMC surface density in normal galaxies
phip=3.0                  # pressure normalized to isothermal gas disk value
vdisp1=0.8e6                # ISM velocity dispersion in local galaxies
vdisp2=3e6                # ISM velocity dispersion in high-z galaxies and SBs
q=1.0                     # Toomre Q
betasb=1.0                # Index of rotation curve in SBs (solid body)
betadisk=0.0              # Index of rotation curve in normal disks (flat)


# Heiderman data
area_heid=[]          # Cloud area in pc^2
sigmag_heid=[]        # Gas surface density
sigmasfr_heid=[]      # SFR surface density
mass_heid=[]          # Total mass
fp=open('heiderman.txt', 'r')
for line in fp:
    spl=line.split()
    area_heid.append(float(spl[8]))
    sigmag_heid.append(float(spl[14]))
    sigmasfr_heid.append(float(spl[20]))
    mass_heid.append(float(spl[11]))
fp.close()
area_heid=array(area_heid)*pc**2
sigmag_heid=array(sigmag_heid)*msun/pc**2
sigmasfr_heid=array(sigmasfr_heid)*msun/pc**2/myr
mass_heid=array(mass_heid)*msun
rho_heid=3*pi**0.5/4*mass_heid/area_heid**1.5
tff_heid=sqrt(3*pi/(32*g*rho_heid))
eff_heid=sigmasfr_heid/(sigmag_heid/tff_heid)


# Lada data; masses for A_K = 0.1 and 0.8 contours
mass_lada1=array([67714., 71828., 99930., 18438., 14964., 14165., 1137., 
                  7937., 2157., 1379., 787.])*msun
mass_lada2=array([13721., 7261., 3199., 1880., 1766., 1296., 258., 178., 
                  163., 124., 75.])*msun
sfr_lada=array([715., 159., 70., 150., 84., 79., 25., 5., 17., 3., 3.])*msun/myr
sigmag_lada2=116.*2.*msun/pc**2
sigmag_lada1=sigmag_lada2/8.0
area_lada1=mass_lada1/sigmag_lada1
area_lada2=mass_lada2/sigmag_lada2
sigmasfr_lada1=sfr_lada/area_lada1
sigmasfr_lada2=sfr_lada/area_lada2/3.0
rho_lada1=3*pi**0.5/4*mass_lada1/area_lada1**1.5
rho_lada2=3*pi**0.5/4*mass_lada2/area_lada2**1.5
tff_lada1=sqrt(3*pi/(32*g*rho_lada1))
tff_lada2=sqrt(3*pi/(32*g*rho_lada2))
eff_lada1=sigmasfr_lada1/(sigmag_lada1/tff_lada1)
eff_lada2=sigmasfr_lada2/(sigmag_lada2/tff_lada2)


# Genzel data
sigmag_genzel=[]
sigmasfr_genzel=[]
sigmagtdyn_genzel=[]
sb_genzel=[]
fp=open('genzel_ks.dat', 'r')
for line in fp:
    if line[0]=='#':
        continue
    spl=line.split()
    if spl[0]=='Name':
        continue
    sigmag_genzel.append(10.0**float(spl[11])*msun/pc**2)
    sigmagtdyn_genzel.append(10.0**float(spl[12])/(2.0*pi)*msun/pc**2/myr)
    sigmasfr_genzel.append(10.0**float(spl[13])*msun/pc**2/myr)
    # Apply CO scalings
    if spl[0][0:3]=='SMM':
        # sub-mm galaxy, so use starburst scalings; Genzel uses alpha = 1.0
        # for starbursts, and Daddi's convention is 0.8
        sb_genzel.append(True)
        sigmag_genzel[-1] = sigmag_genzel[-1]*0.8/1.0*alphasb
        sigmagtdyn_genzel[-1] = sigmagtdyn_genzel[-1]*0.8/1.0*alphasb
    else:
        # non-sub-mm galaxy, so use normal high z scalings; Genzel uses
        # alpha = 3.2, Daddi uses 3.6
        sb_genzel.append(False)
        sigmag_genzel[-1] = sigmag_genzel[-1]*3.6/3.2*alphahiz
        sigmagtdyn_genzel[-1] = sigmagtdyn_genzel[-1]*3.6/3.2*alphahiz
fp.close()
sigmag_genzel=array(sigmag_genzel)
sigmasfr_genzel=array(sigmasfr_genzel)
sigmagtdyn_genzel=array(sigmagtdyn_genzel)
tdyn_genzel=sigmag_genzel/sigmagtdyn_genzel
sb_genzel=array(sb_genzel)
tffgmc_genzel=pi**0.25*vdisp2/(sqrt(8.0)*g*(sigmagmc**3*sigmag_genzel)**0.25)
beta_genzel=betadisk+sigmag_genzel*0
beta_genzel[sb_genzel]=betasb
tfft_genzel=sqrt(3*pi**2*q**2/(32*(beta_genzel+1)*phip))*(tdyn_genzel/(2*pi))
tff_genzel=minimum(tfft_genzel, tffgmc_genzel)
eff_genzel=sigmasfr_genzel/(sigmag_genzel/tff_genzel)


# Bouche
sigmag_bouche=[]
sigmasfr_bouche=[]
sigmagtdyn_bouche=[]
fp=open('KS_2_Bouche.dat', 'r')
for line in fp:
    if line[0]=='#':
        continue
    spl=line.split()
    sigmag_bouche.append(10.0**float(spl[0])*msun/pc**2*alphasb)
    sigmagtdyn_bouche.append(10.0**float(spl[2])/(2.0*pi)*msun/pc**2/myr * 
                             alphasb)
    sigmasfr_bouche.append(10.0**float(spl[1])*msun/pc**2/myr)
fp.close()
sigmag_bouche=array(sigmag_bouche)
sigmasfr_bouche=array(sigmasfr_bouche)
sigmagtdyn_bouche=array(sigmagtdyn_bouche)
tdyn_bouche=sigmag_bouche/sigmagtdyn_bouche
tffgmc_bouche=pi**0.25*vdisp2/(sqrt(8.0)*g*(sigmagmc**3*sigmag_bouche)**0.25)
beta_bouche=betasb
tfft_bouche=sqrt(3*pi**2*q**2/(32*(beta_bouche+1)*phip))*(tdyn_bouche/(2*pi))
tff_bouche=minimum(tfft_bouche, tffgmc_bouche)
eff_bouche=sigmasfr_bouche/(sigmag_bouche/tff_bouche)


# Daddi z = 0.5 data
sigmag_daddiz05=[]
sigmasfr_daddiz05=[]
sigmagtdyn_daddiz05=[]
fp=open('KS_2_Daddiz05.dat', 'r')
for line in fp:
    if line[0]=='#':
        continue
    if line[0]=='\n':
        continue
    spl=line.split()
    sigmag_daddiz05.append(10.0**float(spl[1])*msun/pc**2*alphahiz)
    sigmagtdyn_daddiz05.append(10.0**float(spl[2])*msun/pc**2/myr * 
                               alphahiz)
    sigmasfr_daddiz05.append(10.0**float(spl[0])*msun/pc**2/myr)
fp.close()
sigmag_daddiz05=array(sigmag_daddiz05)
sigmasfr_daddiz05=array(sigmasfr_daddiz05)
sigmagtdyn_daddiz05=array(sigmagtdyn_daddiz05)
tdyn_daddiz05=sigmag_daddiz05/sigmagtdyn_daddiz05
tffgmc_daddiz05=pi**0.25*vdisp2 / \
                 (sqrt(8.0)*g*(sigmagmc**3*sigmag_daddiz05)**0.25)
beta_daddiz05=betadisk
tfft_daddiz05=sqrt(3*pi**2*q**2 / 
                   (32*(beta_daddiz05+1)*phip))*(tdyn_daddiz05/(2*pi))
tff_daddiz05=minimum(tfft_daddiz05, tffgmc_daddiz05)
eff_daddiz05=sigmasfr_daddiz05/(sigmag_daddiz05/tff_daddiz05)


# Daddi z = 2 data
sigmag_daddiz2=[]
sigmasfr_daddiz2=[]
sigmagtdyn_daddiz2=[]
fp=open('KS_2_Daddi.dat', 'r')
for line in fp:
    if line[0]=='#':
        continue
    if line[0]=='\n':
        continue
    spl=line.split()
    sigmag_daddiz2.append(float(spl[0])*msun/pc**2*alphahiz)
    sigmagtdyn_daddiz2.append(float(spl[2])*msun/pc**2/myr * 
                              alphahiz)
    sigmasfr_daddiz2.append(float(spl[1])*msun/pc**2/myr)
fp.close()
sigmag_daddiz2=array(sigmag_daddiz2)
sigmasfr_daddiz2=array(sigmasfr_daddiz2)
sigmagtdyn_daddiz2=array(sigmagtdyn_daddiz2)
tdyn_daddiz2=sigmag_daddiz2/sigmagtdyn_daddiz2
tffgmc_daddiz2=pi**0.25*vdisp2 / \
                 (sqrt(8.0)*g*(sigmagmc**3*sigmag_daddiz2)**0.25)
beta_daddiz2=betadisk
tfft_daddiz2=sqrt(3*pi**2*q**2 / 
                   (32*(beta_daddiz2+1)*phip))*(tdyn_daddiz2/(2*pi))
tff_daddiz2=minimum(tfft_daddiz2, tffgmc_daddiz2)
eff_daddiz2=sigmasfr_daddiz2/(sigmag_daddiz2/tff_daddiz2)


# Tacconi data
sigmag_tacconi=[]
sigmasfr_tacconi=[]
sigmagtdyn_tacconi=[]
fp=open('Tacconi_KS.dat', 'r')
for line in fp:
    if line[0]=='#':
        continue
    if line[0]=='\n':
        continue
    spl=line.split()
    sigmag_tacconi.append(10.0**float(spl[1])*msun/pc**2*alphahiz)
    sigmagtdyn_tacconi.append(10.0**float(spl[2])*msun/pc**2/myr * 
                              alphahiz)
    sigmasfr_tacconi.append(10.0**float(spl[0])*msun/pc**2/myr)
fp.close()
sigmag_tacconi=array(sigmag_tacconi)
sigmasfr_tacconi=array(sigmasfr_tacconi)
sigmagtdyn_tacconi=array(sigmagtdyn_tacconi)
tdyn_tacconi=sigmag_tacconi/sigmagtdyn_tacconi
tffgmc_tacconi=pi**0.25*vdisp2 / \
                 (sqrt(8.0)*g*(sigmagmc**3*sigmag_tacconi)**0.25)
beta_tacconi=betadisk
tfft_tacconi=sqrt(3*pi**2*q**2 / 
                   (32*(beta_tacconi+1)*phip))*(tdyn_tacconi/(2*pi))
tff_tacconi=minimum(tfft_tacconi, tffgmc_tacconi)
eff_tacconi=sigmasfr_tacconi/(sigmag_tacconi/tff_tacconi)


# Kennicutt ULIRG data
sigmag_kenn_ulirg=[]
sigmasfr_kenn_ulirg=[]
tdyn_kenn_ulirg=[]
fp=open('KS_2_KennUlirgs.dat', 'r')
for line in fp:
    if line[0]=='#':
        continue
    if line[0]=='\n':
        continue
    spl=line.split()
    if spl[2]=='-1':
        continue
    sigmag_kenn_ulirg.append(10.0**float(spl[0])*msun/pc**2*alphasb)
    sigmasfr_kenn_ulirg.append(10.0**float(spl[1])*msun/pc**2/myr)
    tdyn_kenn_ulirg.append(float(spl[2])*100*myr)
fp.close()
sigmag_kenn_ulirg=array(sigmag_kenn_ulirg)
sigmasfr_kenn_ulirg=array(sigmasfr_kenn_ulirg)
tdyn_kenn_ulirg=array(tdyn_kenn_ulirg)
tffgmc_kenn_ulirg=pi**0.25*vdisp2 / \
                 (sqrt(8.0)*g*(sigmagmc**3*sigmag_kenn_ulirg)**0.25)
beta_kenn_ulirg=betasb
tfft_kenn_ulirg=sqrt(3*pi**2*q**2 / 
                   (32*(beta_kenn_ulirg+1)*phip))*(tdyn_kenn_ulirg/(2*pi))
tff_kenn_ulirg=minimum(tfft_kenn_ulirg, tffgmc_kenn_ulirg)
eff_kenn_ulirg=sigmasfr_kenn_ulirg/(sigmag_kenn_ulirg/tff_kenn_ulirg)

# Kennicutt spiral data
sigmag_kenn_spiral=[]
sigmasfr_kenn_spiral=[]
tdyn_kenn_spiral=[]
fp=open('KS_2_KennSpirals.dat', 'r')
for line in fp:
    if line[0]=='#':
        continue
    if line[0]=='\n':
        continue
    spl=line.split()
    if spl[2]=='-1':
        continue
    sigmag_kenn_spiral.append(10.0**float(spl[0])*msun/pc**2*alphasb)
    sigmasfr_kenn_spiral.append(10.0**float(spl[1])*msun/pc**2/myr)
    tdyn_kenn_spiral.append(float(spl[2])*100*myr)
fp.close()
sigmag_kenn_spiral=array(sigmag_kenn_spiral)
sigmasfr_kenn_spiral=array(sigmasfr_kenn_spiral)
tdyn_kenn_spiral=array(tdyn_kenn_spiral)
tffgmc_kenn_spiral=pi**0.25*vdisp1 / \
                 (sqrt(8.0)*g*(sigmagmc**3*sigmag_kenn_spiral)**0.25)
beta_kenn_spiral=betadisk
tfft_kenn_spiral=sqrt(3*pi**2*q**2 / 
                   (32*(beta_kenn_spiral+1)*phip))*(tdyn_kenn_spiral/(2*pi))
tff_kenn_spiral=minimum(tfft_kenn_spiral, tffgmc_kenn_spiral)
eff_kenn_spiral=sigmasfr_kenn_spiral/(sigmag_kenn_spiral/tff_kenn_spiral)


# Tacconi 2012 data
name_tacconi12=[]
sigmag_tacconi12=[]
sigmasfr_tacconi12=[]
tdyn_tacconi12=[]
sb_tacconi12=[]
fp=open('tacconi12.txt', 'r')
fp.readline()
fp.readline()
for line in fp:
    spl=line.split()
    if spl[2]=='...':
        continue
    name_tacconi12.append(spl[0])
    sigmag_tacconi12.append(10.**float(spl[-2])*msunpc2)
    sigmasfr_tacconi12.append(10.**float(spl[-1])*msunpc2myr)
    tdyn_tacconi12.append(2*pi*float(spl[3])*kpc/(float(spl[2])*1e5))
    sb_tacconi12.append('merger' in spl[1].lower() or 'amor' in spl[1].lower())
fp.close()
sigmag_tacconi12=array(sigmag_tacconi12)
sigmasfr_tacconi12=array(sigmasfr_tacconi12)
tdyn_tacconi12=array(tdyn_tacconi12)
sb_tacconi12=array(sb_tacconi12)
beta_tacconi12=zeros(len(sb_tacconi12))
beta_tacconi12[sb_tacconi12] = betasb
tffgmc_tacconi12=pi**0.25*vdisp2 / \
                 (sqrt(8.0)*g*(sigmagmc**3*sigmag_tacconi12)**0.25)
tfft_tacconi12=sqrt(3*pi**2*q**2 / 
                   (32*(beta_tacconi12+1)*phip))*(tdyn_tacconi12/(2*pi))
tff_tacconi12=minimum(tfft_tacconi12, tffgmc_tacconi12)
eff_tacconi12=sigmasfr_tacconi12/(sigmag_tacconi12/tff_tacconi12)


# Davis+ 2014 data
data_davis14=ascii.read('DAVIS14_WISEsfrs.txt')
sigmag_davis14=np.array(data_davis14['Sig_gas'])*msunpc2
sigmasfr_davis14=np.array(data_davis14['SFR'])*msunpc2myr
tff_davis14=pi**0.25*vdisp1 / \
             (sqrt(8.0)*g*(sigmagmc**3*sigmag_davis14)**0.25)
eff_davis14=sigmasfr_davis14/(sigmag_davis14/tff_davis14)

# Read in Leroy+ 2013 HERACLES data
heraclesdata=loadtxt('heracles_bins.txt', skiprows=3)
h2_herabin=heraclesdata[:-1,0]
sf_herabin=heraclesdata[:-1,1]
sfscat_herabin=heraclesdata[:-1,2]

# read raw data
fp=open('compile_lit.txt', 'r')
h2_hera=[]
sf_hera=[]
for line in fp:
    if line.strip()[0] == '#':
        continue
    spl=line.split()
    if spl[-1] == 'HERACLES':
        h2_hera.append(float(spl[0]))
        sf_hera.append(float(spl[1]))
h2_hera=np.array(h2_hera)
sf_hera=np.array(sf_hera)
tff_hera=pi**0.25*vdisp1/(sqrt(8.0)*g*(sigmagmc**3*10.**h2_hera*msunpc2)**0.25)
logsigmatff_hera=log10(10.**h2_hera*msunpc2/tff_hera/msunpc2myr)
nbins=100
hist_hera, \
    xedge, \
    yedge = \
              histogram2d(logsigmatff_hera, sf_hera, bins=[nbins,nbins],
                          range=[[-2.0,5.0], [-7,0]])
hist_hera=hist_hera/np.amax(hist_hera)
bbox_hera=[xedge[0], xedge[-1], yedge[0], yedge[-1]]
nrm=colors.Normalize(vmin=-2, vmax=0, clip=True)
img_hera=np.ones((hist_hera.shape[0], hist_hera.shape[1], 3))
img_hera[:,:,0]=1-nrm(np.transpose(np.log10(hist_hera+1e-10)))
img_hera[:,:,1]=1-nrm(np.transpose(np.log10(hist_hera+1e-10)))
#img_hera=np.zeros((hist_hera.shape[0], hist_hera.shape[1], 3))
#img_hera[:,:,2]=nrm(np.transpose(np.log10(hist_hera+1e-10)))

# Bolatto data
SigmaSFRH2_bol_12=np.loadtxt('fig2_12pc.dat', skiprows=2)
SigmaSFRH2_bol_12[np.isnan(SigmaSFRH2_bol_12)] = 0.0
SigmaSFRH2_bol_interp = np.zeros((nbins,nbins))
fp=open('fig2_12pc.dat')
tmp=fp.readline().split()
tmp1=[float(t) for t in tmp]
logSigmaH2_bol_12=np.array(tmp1)
tmp=fp.readline().split()
tmp1=[float(t) for t in tmp]
logSigmaSFR_bol_12=np.array(tmp1)
fp.close()
tff_bol_12=pi**0.25*vdisp1 / \
            (sqrt(8.0)*g*(sigmagmc**3*10.**logSigmaH2_bol_12*msunpc2)**0.25)
logsigmatff_bol=log10(10.**logSigmaH2_bol_12*msunpc2/tff_bol_12/msunpc2myr)
bbox=(logsigmatff_bol[0]-0.5*(logsigmatff_bol[1]-logsigmatff_bol[0]),
      logsigmatff_bol[-1]+0.5*(logsigmatff_bol[1]-logsigmatff_bol[0]),
      logSigmaSFR_bol_12[0]-0.5*(logSigmaSFR_bol_12[1]-logSigmaSFR_bol_12[0]),
      logSigmaSFR_bol_12[-1]+0.5*(logSigmaSFR_bol_12[1]-logSigmaSFR_bol_12[0]))
interpshift(SigmaSFRH2_bol_12, bbox, SigmaSFRH2_bol_interp, bbox_hera)
SigmaSFRH2_bol_interp=SigmaSFRH2_bol_interp/np.amax(SigmaSFRH2_bol_interp)
img_hera[:,:,0]=img_hera[:,:,0]-nrm(log10(SigmaSFRH2_bol_interp+1e-10))
img_hera[:,:,2]=img_hera[:,:,2]-nrm(log10(SigmaSFRH2_bol_interp+1e-10))
#img_hera[:,:,1]=nrm(log10(SigmaSFRH2_bol_interp+1e-10))

# Suppress looping around in color space
img_hera[:,:,0]=img_hera[:,:,0]*(img_hera[:,:,0]>=0)


# Evans 2013 data
evansdat=ascii.read('logsdens.dat')
logsigmagtff_evans=evansdat['col1']
logsigmasfr_evans=evansdat['col4']
evansdat1=ascii.read('logsdens1.dat')
logsigmagtff_evans1=evansdat1['col1']
logsigmasfr_evans1=evansdat1['col4']
logsigmasfrerr_evans1=evansdat1['col6']
evansdat2=ascii.read('logsdensu.dat')
logsigmagtff_evans2=evansdat2['col1']
logsigmasfr_evans2=evansdat2['col4']
logsigmasfrerr_evans2=evansdat2['col6']

# Lada 2013 / Lombardi data
ladadat=ascii.read('lombardi3.dat')
arealada=ladadat['AREA'].data*pc**2
sigmaSFRlada=ladadat['NSTARS'].data*0.25*msun/(0.25*myr)/arealada
sigmaSFRerrlada=sqrt(ladadat['NSTARS'].data)*0.25*msun/(0.25*myr)/arealada
AKtosigmag = 183*msun/pc**2
sigmaglada=ladadat['AK_MEAN']*AKtosigmag
masslada=sigmaglada*arealada
areadifflada=arealada[:-1]-arealada[1:]
massdifflada=masslada[:-1]-masslada[1:]
sigmaSFRdifflada=(ladadat['NSTARS'].data[:-1] - 
                  ladadat['NSTARS'].data[1:])*0.25*msun/(0.25*myr) \
    / areadifflada
sigmaSFRdifferrlada=sqrt(ladadat['NSTARS'].data[:-1] - 
                         ladadat['NSTARS'].data[1:]) \
    *0.25*msun/(0.25*myr) \
    / areadifflada
sigmagdifflada=massdifflada/areadifflada
#rholada=3*pi**0.5/4*masslada/arealada**1.5
rholada=3*pi**0.5/4*massdifflada/areadifflada**1.5
tfflada=sqrt(3*pi/(32*g*rholada))
#efflada=sigmaSFRlada/(sigmaglada/tfflada)
efflada=sigmaSFRdifflada/(sigmagdifflada/tfflada)


# Gutermuth 2011 data
gutermuthdata=ascii.read('gutermuth11.txt')
area_gutermuth11=gutermuthdata['Area[pc^2]']*pc**2
mass_gutermuth11=gutermuthdata['Mass[msol]']*msun
sigmag_gutermuth11=mass_gutermuth11/area_gutermuth11
nstar_gutermuth11=gutermuthdata['N_II']
sigmasfr_gutermuth11=nstar_gutermuth11*0.5*msun/(2.0*myr)/area_gutermuth11
rho_gutermuth11 = 3.*sqrt(pi)/4.*sigmag_gutermuth11 * \
                   area_gutermuth11**(-0.5)
tff_gutermuth11 = sqrt(3.*pi/(32.*g*rho_gutermuth11))

# Wu 2010 data
wudata=ascii.read('wu10.txt')
area_wu10=wudata['Area']*pc**2
sigmag_wu10=10.**wudata['logSig_gas']*msunpc2
sigmasfr_wu10=10.**wudata['logSig_SFR']*msunpc2myr
rho_wu10 = 3.*sqrt(pi)/4.*sigmag_wu10 * \
                   area_wu10**(-0.5)
tff_wu10 = sqrt(3.*pi/(32.*g*rho_wu10))

# Make plot

# Set general paramters
ms=8
mew=2

# Open figure window
fig=plt.figure(1, figsize=(8,8))
ax=plt.subplot(111)
plt.subplots_adjust(left=0.125, right=0.8, bottom=0.1, top=0.9)

# Heracles and Bolatto data
imshow(img_hera, interpolation='nearest', aspect='auto', origin='lower',
       extent=bbox_hera)

# e_ff = 0.01
p1,=plot(arange(-2,5.1,0.1), arange(-2,5.1,0.1)+log10(0.01), 'k', lw=5)
fill_between(arange(-2,5.1,0.1), arange(-2,5.1,0.1)+log10(0.01*3),
             arange(-2,5.1,0.1)+log10(0.01/3), color='k', 
             alpha=0.25, lw=0)

# Milky Way data
p_heid,=plot(log10(sigmag_heid/tff_heid/msunpc2myr), 
             log10(sigmasfr_heid/msunpc2myr), 's', 
             mec='k', ms=ms, color='#ffaaaa')
p_lada,=plot(log10(sigmag_lada1/tff_lada1/msunpc2myr), 
             log10(sigmasfr_lada1/msunpc2myr), 'o', 
             mec='k', ms=ms, color='#ff7777')
plot(log10(sigmag_lada2/tff_lada2/msunpc2myr), 
     log10(sigmasfr_lada2/msunpc2myr), 'o', 
     mec='k', ms=ms, color='#ff7777')
yerr=zeros((2,sigmag_wu10.shape[0]))
yerr[1,:]=0.75
errorbar(log10(sigmag_wu10/tff_wu10/msunpc2myr),
         log10(sigmasfr_wu10/msunpc2myr),
         yerr=yerr, uplims=True,
         color='#ff0000', ms=ms, fmt='*', ecolor='#ff0000', lw=2,
         capthick=2)
p_wu = Line2D(range(1), range(1), color="w", marker='*', mfc="#ff0000", ms=ms)
p_gutermuth,=plot(log10(sigmag_gutermuth11/tff_gutermuth11/msunpc2myr),
     log10(sigmasfr_gutermuth11/msunpc2myr),
                  '^', color='#dd0000', ms=ms, mec='k')
p_lada1,=plot(np.log10(sigmagdifflada/tfflada/msunpc2myr), 
              np.log10(sigmaSFRdifflada/msunpc2myr), 'D', 
              color='#bb0000', mec='k', 
              ms=ms)
p_evans1,=plot(logsigmagtff_evans, logsigmasfr_evans, 'p', 
               color='#aa0000', mec='k', ms=ms)
#plot(logsigmagtff_evans1, logsigmasfr_evans1+logsigmasfrerr_evans1, 
#     'rv', mec='k', ms=ms)
#plot(logsigmagtff_evans2, logsigmasfr_evans2+logsigmasfrerr_evans2, 
#     'rv', mec='k', ms=ms)
yerr=zeros((2,logsigmasfrerr_evans1.shape[0]))
yerr[0,:]=2*logsigmasfrerr_evans1
errorbar(logsigmagtff_evans1, logsigmasfr_evans1+logsigmasfrerr_evans1, 
         yerr=yerr, lolims=True,
         color='#aa0000', ms=ms, fmt='p', ecolor='#aa0000', lw=2,
         capthick=2)
yerr=zeros((2,logsigmasfrerr_evans2.shape[0]))
yerr[0,:]=2*logsigmasfrerr_evans2
errorbar(logsigmagtff_evans2, logsigmasfr_evans2+logsigmasfrerr_evans2, 
         yerr=yerr, lolims=True,
         color='#aa0000', ms=ms, fmt='p', ecolor='#aa0000', lw=2,
         capthick=2)

# z = 0 whole-galaxy data
#plot(log10(sigmag_kenn_ulirg/tff_kenn_ulirg/msunpc2myr), 
#     log10(sigmasfr_kenn_ulirg/msunpc2myr), 
#     'gs', mfc='none', mec='g', ms=ms, mew=mew)
plot(log10(sigmag_kenn_ulirg/tff_kenn_ulirg/msunpc2myr), 
     log10(sigmasfr_kenn_ulirg/msunpc2myr), 
     'gs', mec='k', ms=ms, mfc='#00d000')
p_kenn,=plot(log10(sigmag_kenn_spiral/tff_kenn_spiral/msunpc2myr), 
     log10(sigmasfr_kenn_spiral/msunpc2myr), 
     'gs', ms=ms, mec='k')
p_davis,=plot(log10(sigmag_davis14/tff_davis14/msunpc2myr),
              log10(sigmasfr_davis14/msunpc2myr),
              'o', mfc='#00a000', ms=ms, mec='k')

# High-z data
highzcolor='m'
#plot(log10(sigmag_tacconi/tff_tacconi/msunpc2myr), 
#     log10(sigmasfr_tacconi/msunpc2myr), highzcolor+'*', 
#     ms=ms+1, mec=highzcolor)
p_tacc,=plot(log10(sigmag_tacconi12/tff_tacconi12/msunpc2myr), 
     log10(sigmasfr_tacconi12/msunpc2myr), highzcolor+'*', 
     ms=ms+1, mec='k')
p_daddi,=plot(log10(sigmag_daddiz05/tff_daddiz05/msunpc2myr), 
     log10(sigmasfr_daddiz05/msunpc2myr), 
     highzcolor+'o', ms=ms, mec='k')
plot(log10(sigmag_daddiz2/tff_daddiz2/msunpc2myr), 
     log10(sigmasfr_daddiz2/msunpc2myr), 
     highzcolor+'o', ms=ms, mec='k')
#plot(log10(sigmag_genzel[sb_genzel]/tff_genzel[sb_genzel]/msunpc2myr),
#     log10(sigmasfr_genzel[sb_genzel]/msunpc2myr),
#     highzcolor+'p', mfc='none', mew=mew, ms=ms, mec=highzcolor)
plot(log10(sigmag_genzel[sb_genzel]/tff_genzel[sb_genzel]/msunpc2myr),
     log10(sigmasfr_genzel[sb_genzel]/msunpc2myr),
     highzcolor+'p', mfc='#ff88ff', ms=ms, mec='k')
p_genzel,=plot(log10(sigmag_genzel[logical_not(sb_genzel)] / 
     tff_genzel[logical_not(sb_genzel)]/msunpc2myr),
     log10(sigmasfr_genzel[logical_not(sb_genzel)]/msunpc2myr),
     highzcolor+'p', ms=ms, mec='k')
#plot(log10(sigmag_bouche/tff_bouche/msunpc2myr), 
#     log10(sigmasfr_bouche/msunpc2myr), 
#     highzcolor+'s', ms=ms, mec=highzcolor, mfc='none', mew=mew)
p_bouche,=plot(log10(sigmag_bouche/tff_bouche/msunpc2myr), 
     log10(sigmasfr_bouche/msunpc2myr), 
     highzcolor+'s', ms=ms, mec='k', mfc='#ff88ff')

# Set range
xlim([-2,5])
ylim([-5,4])

# Add labels
xlabel(r'$\log\,\Sigma_{\rm H_2}/t_{\rm ff}$ $[M_\odot\,{\rm pc}^{-2}\,{\rm Myr}^{-1}]$')
ylabel(r'$\log\,\Sigma_{\rm SFR}$ $[M_\odot\,{\rm pc}^{-2}\,{\rm Myr}^{-1}]$')

# Legend
prd = Rectangle((0,0), 1, 1, fc='r')
pgr = Rectangle((0,0), 1, 1, fc='#00ff00')
pbl = Rectangle((0,0), 1, 1, fc='b')
blank = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', 
                  linewidth=0)
l1 = legend([p_heid, p_lada, p_wu, p_gutermuth, p_lada1, p_evans1], 
            ["Heiderman+ 2010",
             "Lada+ 2010",
             "Wu+ 2010",
             "Gutermuth+ 2011",
             "Lada+ 2013",
             "Evans+ 2014"],
            loc=2, numpoints=1,
            title="Milky Way")
l2 = legend([p_kenn, p_davis, p_bouche, p_daddi, p_genzel, p_tacc],
            [
             "Kennicutt 1998",
             "Davis+ 2014",
             "Bouche+ 2007",
             "Daddi+ 2008, 2010",
             "Genzel+ 2010",
             "Tacconi+ 2013"],
            loc=4, numpoints=1, title='Unresolved galaxies')
l3 = legend([p1], [r"$\epsilon_{\rm ff} = 0.01$"], loc=1)
l4 = legend([pbl, pgr], ["Leroy+ 2013",
                         "Bolatto+ 2011"], loc='lower left',
            title='Resolved galaxies')
ax.add_artist(l1)
ax.add_artist(l2)
ax.add_artist(l3)

# Colorbar setup
cdictbl={ 'red' : [ (0.0, 1.0, 1.0),
                    (1.0, 0.0, 0.0) ],
          'green' : [ (0.0, 1.0, 1.0),
                      (1.0, 0.0, 0.0) ],
          'blue' : [ (0.0, 1.0, 1.0),
                     (1.0, 1.0, 1.0) ] }
cdictgr={ 'red' : [ (0.0, 1.0, 1.0),
                    (1.0, 0.0, 0.0) ],
          'green' : [ (0.0, 1.0, 1.0),
                      (1.0, 1.0, 1.0) ],
          'blue' : [ (0.0, 1.0, 1.0),
                     (1.0, 0.0, 0.0) ] }
cdictrd={ 'red' : [ (0.0, 1.0, 1.0),
                    (1.0, 1.0, 1.0) ],
          'green' : [ (0.0, 1.0, 1.0),
                      (1.0, 0.0, 0.0) ],
          'blue' : [ (0.0, 1.0, 1.0),
                     (1.0, 0.0, 0.0) ] }
cmapbl=colors.LinearSegmentedColormap('bl', cdictbl)
cmapgr=colors.LinearSegmentedColormap('gr', cdictgr)
cmaprd=colors.LinearSegmentedColormap('rd', cdictrd)
cbarwidth=0.025

# Colorbars
pts=ax.get_axes().get_position().get_points()
axcbar1=fig.add_axes([pts[1,0], pts[0,1], cbarwidth, pts[1,1]-pts[0,1]], label='bar1')
norm=colors.Normalize(vmin=-2, vmax=0.0)
cbar.ColorbarBase(axcbar1, norm=norm, orientation='vertical', cmap=cmapbl)
plt.setp(axcbar1.get_yticklabels(), visible=False)
axcbar2=fig.add_axes([pts[1,0]+cbarwidth, pts[0,1], cbarwidth, pts[1,1]-pts[0,1]], label='bar2')
cbar.ColorbarBase(axcbar2, norm=norm, orientation='vertical', 
                  ticks=[-2,-1.5,-1.0,-0.5,0], cmap=cmapgr)
axcbar2.set_ylabel('Log point density')
axcbar2.yaxis.set_ticks_position('right')
axcbar2.yaxis.set_label_position('right')

# Save
#savefig("effplot.pdf")
