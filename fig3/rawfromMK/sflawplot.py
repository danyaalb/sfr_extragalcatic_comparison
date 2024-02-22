try:
    import astropy.io.ascii as asciitable
except ImportError:
    import asciitable
import numpy as np
import numpy.ma as ma
import collections
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.colorbar as cbar
import scipy.interpolate as interp
from interpshift import interpshift

Msun=1.99e33
pc=3.09e18
yr=365.25*24.*3600.
Myr=1e6*yr
G=6.67e-8
sigmaISM=8e5

nbins=100

# Bigiel data
try:
    indata
except:
    indata=asciitable.read('bigiel1.txt', Reader=asciitable.Cds)
    outdata=np.loadtxt('things_sflaw_outer_krumholz.txt', delimiter='&',
                       skiprows=18, usecols=(1,2))

    # Count number of points per galaxy to set weight
    fp=open('things_sflaw_outer_krumholz.txt', 'r')
    for i in range(18):
        fp.readline()
    galnameout=[]
    for line in fp:
        galnameout.append(line.partition('&')[0].strip())
    fp.close()
    outnamectr=collections.Counter(galnameout)
    outwgt=np.zeros(len(outdata))
    for n in outnamectr:
        i1=galnameout.index(n)
        i2=len(galnameout)-galnameout[::-1].index(n)
        outwgt[i1:i2] = 1.0/outnamectr[n]
    galnamein=indata['Name'].tolist()
    innamectr=collections.Counter(galnamein)
    inwgt=np.zeros(len(indata))
    for n in innamectr:
        i1=galnamein.index(n)
        i2=len(galnamein)-galnamein[::-1].index(n)
        inwgt[i1:i2] = 1.0/innamectr[n]

    # Get histogram of inner galaxy data
    hiin=ma.filled(10.**indata['logHI'], 1e-50)
    h2in=ma.filled(10.**indata['logH2'], 1e-50)
    gasin=hiin+h2in
    loggasin=np.log10(gasin)
    logh2in=np.log10(h2in)
    logsfrin=ma.filled(indata['logSFR'], -50)
    histin, xedgein, \
        yedgein = \
        np.histogram2d(loggasin, logsfrin,
                       bins=[nbins,nbins],
                       range=[[-0.5,3.1], [-7,0]])
    xcenin=0.5*(xedgein[:-1]+xedgein[1:])
    ycenin=0.5*(yedgein[:-1]+yedgein[1:])
    histin=histin/np.amax(histin)
    histh2in, dummy, \
        dummy = \
        np.histogram2d(logh2in, logsfrin, bins=[nbins,nbins],
                       range=[[-0.5,3.1], [-7,0]])
    histh2in=histh2in/np.amax(histh2in)
    bbox_in=[xedgein[0], xedgein[-1], yedgein[0], yedgein[-1]]

    # Get histogram of outer galaxy data
    loggasout=outdata[:,0]
    logsfrout=outdata[:,1]
    histout, xedgeout, \
        yedgeout = \
        np.histogram2d(loggasout, logsfrout,
                       bins=[nbins,nbins], weights=outwgt,
                       range=[[-0.5,3.1], [-7,0]])
    xcenout=0.5*(xedgeout[:-1]+xedgeout[1:])
    ycenout=0.5*(yedgeout[:-1]+yedgeout[1:])
    histout=histout/np.amax(histout)
    bbox_out=[xedgeout[0], xedgeout[-1], yedgeout[0], yedgeout[-1]]

    # Bigiel's median and scatter
    binctr=np.array([-0.150515, 0.150515, 0.451545, 0.752575])
    binmed=np.array([-5.50282, -5.05890, -4.47417, -3.98886])
    binscatter=np.array([1.5, 1.2, 0.64, 0.46])

# Now read Bolatto et al data
try:
    SigmaSFR_bol_interp
except:
    # Read original Bolatto data
    SigmaSFR_bol=np.loadtxt('logSsfrlogSgas_200pc.dat')
    logSigma_bol=np.loadtxt('logSsfrlogSgas_200pc_x.dat')
    logSFR_bol=np.loadtxt('logSsfrlogSgas_200pc_y.dat')

    # Construct a bivariate spline interpolation over it
    bbox=(logSigma_bol[0]-0.5*(logSigma_bol[1]-logSigma_bol[0]),
            logSigma_bol[-1]+0.5*(logSigma_bol[1]-logSigma_bol[0]),
            logSFR_bol[0]-0.5*(logSFR_bol[1]-logSFR_bol[0]),
            logSFR_bol[-1]+0.5*(logSFR_bol[1]-logSFR_bol[0]))
    bbox_1=(logSigma_bol[0], logSigma_bol[-1], logSFR_bol[0], logSFR_bol[-1])
    bolinterp=interp.RectBivariateSpline(logSigma_bol, logSFR_bol, 
                                         SigmaSFR_bol, bbox=bbox,
                                         kx=1, ky=1)

    # Within the range of Bolatto's data, interpolate onto Bigiel's grid
    xidx=np.where(np.logical_and(xcenin > bbox[0], xcenin < bbox[1]))[0]
    yidx=np.where(np.logical_and(ycenin > bbox[2], ycenin < bbox[3]))[0]
    xcoord=xcenin[xidx]
    ycoord=ycenin[yidx]
    SigmaSFR_bol_interp = np.zeros((nbins,nbins))
    #SigmaSFR_bol_interp[xidx[0]-1:xidx[-1], yidx[0]-1:yidx[-1]] = \
    #    bolinterp(xcoord, ycoord)
    bbox_new=(xcenin[0]-0.5*(xcenin[1]-xcenin[0]),
              xcenin[-1]+0.5*(xcenin[1]-xcenin[0]),
              ycenin[0]-0.5*(ycenin[1]-ycenin[0]),
              ycenin[-1]+0.5*(ycenin[1]-ycenin[0]))
    interpshift(SigmaSFR_bol, bbox, SigmaSFR_bol_interp, bbox_new)
    bbox_new_1=(xcenin[0], xcenin[-1], ycenin[0], ycenin[-1])
    SigmaSFR_bol_interp=SigmaSFR_bol_interp/np.amax(SigmaSFR_bol_interp)

# Bolatto et al H2 data
try:
    SigmaH2SFR_bol_200
except:
    SigmaH2SFR_bol_200=np.loadtxt('fig2_200pc.dat')
    SigmaH2SFR_bol_1000=np.loadtxt('fig2_kpc.dat')
    fp=open('fig2_12pc.dat')
    tmp=fp.readline().split()
    tmp1=[float(t) for t in tmp]
    logSigmaH2_bol_12=np.array(tmp1)
    tmp=fp.readline().split()
    tmp1=[float(t) for t in tmp]
    logSigmaSFR_bol_12=np.array(tmp1)
    fp.close()
    SigmaSFRH2_bol_12=np.loadtxt('fig2_12pc.dat', skiprows=2)
    SigmaSFRH2_bol_12[np.isnan(SigmaSFRH2_bol_12)] = 0.0
    SigmaSFRH2_bol_interp = np.zeros((nbins,nbins))
    bbox=(logSigmaH2_bol_12[0]-0.5*(logSigmaH2_bol_12[1]-logSigmaH2_bol_12[0]),
            logSigmaH2_bol_12[-1]+0.5*(logSigmaH2_bol_12[1]-logSigmaH2_bol_12[0]),
            logSigmaSFR_bol_12[0]-0.5*(logSigmaSFR_bol_12[1]-logSigmaSFR_bol_12[0]),
            logSigmaSFR_bol_12[-1]+0.5*(logSigmaSFR_bol_12[1]-logSigmaSFR_bol_12[0]))
    interpshift(SigmaSFRH2_bol_12, bbox, SigmaSFRH2_bol_interp, bbox_new)
    SigmaSFRH2_bol_interp=SigmaSFRH2_bol_interp/np.amax(SigmaSFRH2_bol_interp)


# Schruba's data
try:
    sdata
except:
    sdata=np.loadtxt('schruba11.txt', usecols=(1,2,3,4))
    shi=sdata[:,0]
    sh2=sdata[:,1]
    ssfr=sdata[:,2]
    squal=sdata[:,3]

# Rafelski's LGB outskirts and DLA data
SigmaLBG=np.array([86.891123, 81.115226, 67.540202, 70.898418, 52.409843, 
                   57.544882, 51.858513, 40.570607, 51.944577])
SigmaSFRLBG=np.array([0.0064785162, 0.0045893155, 0.0039155767, 0.0023390654,
                      0.0017876233, 0.0016736985, 0.0011322976, 0.0018736368,
                      0.0013240842])
SigmaLBGlow=-np.array([-2.5709212, -1.5190513, -1.1447733, -2.1883832,
                      -2.8490171, -3.8634207, -7.2130007, -1.4388164,
                      -1.8421886])
SigmaLBGhi=-np.array([-2.4004563, -2.8464341, -1.1000290, -2.0374478,
                     -2.5201820, -3.3272074, -5.4068872, -1.3259522,
                     -1.6976829])
SigmaLBGlowlog=-np.log10((-SigmaLBGlow+SigmaLBG)/SigmaLBG)
SigmaLBGhilog=np.log10((SigmaLBGhi+SigmaLBG)/SigmaLBG)
SigmaLBGerr=np.transpose(np.array(zip(SigmaLBGlowlog, SigmaLBGhilog)))
SigmaSFRLBGhilo=np.array([0.00096271353, 0.00080198775, 0.00064601040,
                          0.00069254802, 0.00056371335, 0.00061655062,
                          0.00061522698, 0.00049356358, 0.00055216044])
SigmaSFRLBGlowlog=-np.log10((-SigmaSFRLBGhilo+SigmaSFRLBG)/SigmaSFRLBG)
SigmaSFRLBGhilog=np.log10((SigmaSFRLBGhilo+SigmaSFRLBG)/SigmaSFRLBG)
SigmaSFRLBGerr=np.transpose(np.array(zip(SigmaSFRLBGlowlog, SigmaSFRLBGhilog)))
dataDLA=np.array([[68.726716,     0.022300453],
                  [68.206464,     0.020495037],
                  [70.040084,     0.013797099],
                  [72.472339,     0.011055385],
                  [70.801623,    0.0094356587],
                  [72.157125,    0.0063930909],
                  [59.657986,    0.0048218819],
                   [45.846304,    0.0062466624]])
SigmaDLA=dataDLA[:,0]
SigmaSFRDLA=dataDLA[:,1]


# Leroy '13 HERACLES data
try:
    histh2in_hera
except:
    # read binned data
    heraclesdata=np.loadtxt('heracles_bins.txt', skiprows=3)
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
    histh2in_hera, dummy, \
        dummy = \
        np.histogram2d(h2_hera, sf_hera, bins=[nbins,nbins],
                       range=[[-0.5,3.1], [-7,0]])
    histh2in_hera=histh2in_hera/np.amax(histh2in_hera)


# Construct image
#img=np.zeros((histin.shape[0], histin.shape[1], 3))
img=np.ones((histin.shape[0], histin.shape[1], 3))

# Log scale
nrm=colors.Normalize(vmin=-2, vmax=0, clip=True)
scalefac=sqrt(2)/2
#img[:,:,2]=nrm(np.transpose(np.log10(histin)))
#img[:,:,1]=nrm(np.log10(SigmaSFR_bol_interp))
#img[:,:,0]=nrm(np.transpose(np.log10(histout)))
img[:,:,1]=img[:,:,1]-nrm(np.transpose(np.log10(histin)))*scalefac
img[:,:,0]=img[:,:,0]-nrm(np.transpose(np.log10(histin)))*scalefac
img[:,:,0]=img[:,:,0]-nrm(np.log10(SigmaSFR_bol_interp))*scalefac
img[:,:,2]=img[:,:,2]-nrm(np.log10(SigmaSFR_bol_interp))*scalefac
img[:,:,1]=img[:,:,1]-nrm(np.transpose(np.log10(histout)))*scalefac
img[:,:,2]=img[:,:,2]-nrm(np.transpose(np.log10(histout)))*scalefac
img[img < 0]=0.0

# Linear scale
#img[:,:,2]=np.transpose(np.clip(histin,0,1))
#img[:,:,1]=np.clip(SigmaSFR_bol_interp, 0, 1)
#img[:,:,0]=np.transpose(np.clip(histout,0,1))


# Plot
fig=plt.figure(1, figsize=(8,6))
ax=plt.subplot(111)
plt.subplots_adjust(left=0.125, right=0.8, bottom=0.1, top=0.9)
plt.imshow(img, extent=bbox_in,
           interpolation='nearest', aspect='auto', origin='lower')
p1=plt.errorbar(binctr, binmed, yerr=binscatter, marker='o', ms=15, mfc='r', 
                mec='k', ecolor='r', elinewidth=3, capsize=5, capthick=3, 
                fmt='o')
p2,=plt.plot(np.log10(SigmaLBG), np.log10(SigmaSFRLBG), 'mo', ms=10)
p3,=plt.plot(np.log10(SigmaDLA), np.log10(SigmaSFRDLA), 'mv', ms=10)

# Lines of constant depletion time
xvec=10.**np.arange(-0.5,2.501,0.1)
plt.plot(np.log10(xvec), np.log10(xvec/1e3), color="0.5", linestyle="--",
         linewidth=2)
plt.text(2.45, np.log10(10.**2.45/1e3)-0.17, "1 Gyr", color="0.5", va="bottom", 
         rotation="25", ha="right")
plt.plot(np.log10(xvec), np.log10(xvec/1e4), color="0.5", linestyle="--",
         linewidth=2)
plt.text(2.45, np.log10(10.**2.45/1e4)-0.17, "10 Gyr", color="0.5", va="bottom", 
         rotation="25", ha="right")
plt.plot(np.log10(xvec), np.log10(xvec/1e5), color="0.5", linestyle="--",
         linewidth=2)
plt.text(2.45, np.log10(10.**2.45/1e5)-0.25, "100 Gyr", color="0.5", va="bottom", 
         rotation="25", ha="right")
plt.plot(np.log10(xvec), np.log10(xvec/1e6), color="0.5", linestyle="--",
         linewidth=2)
plt.text(2.45, np.log10(10.**2.45/1e6)-0.3, "1000 Gyr", color="0.5", va="bottom", 
         rotation="25", ha="right")

# Limits and labels
xlim([-0.5,2.5])
ylim([-6,0])
xlabel(r'$\log\,\Sigma$ $[M_\odot\,{\rm pc}^{-2}]$')
ylabel(r'$\log\,\Sigma_{\rm SFR}$ $[M_\odot\,{\rm pc}^{-2}\,{\rm Myr}^{-1}]$')

# Legend
prd = Rectangle((0,0), 1, 1, fc='r')
pgr = Rectangle((0,0), 1, 1, fc='#00ff00')
pbl = Rectangle((0,0), 1, 1, fc='b')
legend([pbl, pgr, prd, p1, p2, p3], 
       ["Inner disks",
        "SMC",
        "Outer disks",
        "Outer disks (median)",
        "LBG outskirts",
        "DLAs (upper limits)"],
       loc=2, numpoints=1)

# Colorbars
cdictbl={ 'red' : [ (0.0, 1.0, 1.0),
                    (1.0, 0.0, 1.0-scalefac) ],
          'green' : [ (0.0, 1.0, 1.0),
                      (1.0, 0.0, 1.0-scalefac) ],
          'blue' : [ (0.0, 0.0, 1.0),
                     (1.0, 1.0, 1.0) ] }
cdictgr={ 'red' : [ (0.0, 0.0, 1.0),
                    (1.0, 0.0, 1.0-scalefac) ],
          'green' : [ (0.0, 0.0, 1.0),
                      (1.0, 1.0, 1.0) ],
          'blue' : [ (0.0, 0.0, 1.0),
                     (1.0, 0.0, 1.0-scalefac) ] }
cdictrd={ 'red' : [ (0.0, 0.0, 1.0),
                    (1.0, 1.0, 1.0) ],
          'green' : [ (0.0, 0.0, 1.0),
                      (1.0, 0.0, 1.0-scalefac) ],
          'blue' : [ (0.0, 0.0, 1.0),
                     (1.0, 0.0, 1.0-scalefac) ] }
cmapbl=colors.LinearSegmentedColormap('bl', cdictbl)
cmapgr=colors.LinearSegmentedColormap('gr', cdictgr)
cmaprd=colors.LinearSegmentedColormap('rd', cdictrd)
cbarwidth=0.025
pts=ax.get_axes().get_position().get_points()
axcbar1=fig.add_axes([pts[1,0], pts[0,1], cbarwidth, pts[1,1]-pts[0,1]], label='bar1')
norm=colors.Normalize(vmin=-2, vmax=0.0)
cbar.ColorbarBase(axcbar1, norm=norm, orientation='vertical', cmap=cmapbl)
plt.setp(axcbar1.get_yticklabels(), visible=False)
axcbar2=fig.add_axes([pts[1,0]+cbarwidth, pts[0,1], cbarwidth, pts[1,1]-pts[0,1]], label='bar2')
cbar.ColorbarBase(axcbar2, norm=norm, orientation='vertical', 
                  ticks=[-2,-1.5,-1.0,-0.5,0], cmap=cmapgr)
plt.setp(axcbar2.get_yticklabels(), visible=False)
axcbar3=fig.add_axes([pts[1,0]+2*cbarwidth, pts[0,1], cbarwidth, pts[1,1]-pts[0,1]], label='bar3')
cbar.ColorbarBase(axcbar3, norm=norm, orientation='vertical', 
                  ticks=[-2,-1.5,-1.0,-0.5,0], cmap=cmaprd)
axcbar3.set_ylabel('Log point density')
axcbar3.yaxis.set_ticks_position('right')
axcbar3.yaxis.set_label_position('right')

#savefig('sflawplot.pdf')


# Next figure: H2 only

# Construct H2 image
#imgh2=np.zeros((histin.shape[0], histin.shape[1], 3))
imgh2=np.ones((histin.shape[0], histin.shape[1], 3))

# Log scale
#imgh2[:,:,2]=nrm(np.transpose(np.log10(histh2in)))
#imgh2[:,:,2]=nrm(np.transpose(np.log10(histh2in_hera)))
#imgh2[:,:,1]=nrm(np.log10(SigmaSFRH2_bol_interp))
imgh2[:,:,0]=imgh2[:,:,0]-nrm(np.transpose(np.log10(histh2in_hera)))*scalefac
imgh2[:,:,1]=imgh2[:,:,1]-nrm(np.transpose(np.log10(histh2in_hera)))*scalefac
imgh2[:,:,0]=imgh2[:,:,0]-nrm(np.log10(SigmaSFRH2_bol_interp))
imgh2[:,:,2]=imgh2[:,:,2]-nrm(np.log10(SigmaSFRH2_bol_interp))
imgh2[imgh2 < 0]=0.0

# Plot
fig=plt.figure(2, figsize=(8,6))
ax=plt.subplot(111)
plt.subplots_adjust(left=0.125, right=0.8, bottom=0.1, top=0.9)
plt.imshow(imgh2, extent=bbox_in,
           interpolation='nearest', aspect='auto', origin='lower')

# Add Schruba data
p1,=plt.plot(np.log10(sh2[squal==2]), np.log10(ssfr[squal==2]), 'ro', ms=6)
p2,=plt.plot(np.log10(sh2[squal==1]), np.log10(ssfr[squal==1]), 'ro', ms=4)

# Add Bolatto data
p3,=plt.plot(SigmaH2SFR_bol_200[0,:], SigmaH2SFR_bol_200[1,:], marker='s', mfc='#00ff00', 
             ms=10, mec='k',
             linewidth=0)
p4,=plt.plot(SigmaH2SFR_bol_1000[:,0], SigmaH2SFR_bol_1000[:,1], marker='o', mfc='#00ff00',
             ms=10, mec='k',
             linewidth=0)

# Heracles bins
p5=plt.errorbar(h2_herabin, sf_herabin, yerr=sfscat_herabin, marker='o', 
                ms=15, mfc='b', 
                mec='k', ecolor='b', elinewidth=3, capsize=5, capthick=3, 
                fmt='o')

# Lines of constant depletion time
xvec=10.**np.arange(-0.5,2.501,0.1)
plt.plot(np.log10(xvec), np.log10(xvec/1e3), color="0.5", linestyle="--",
         linewidth=2)
plt.text(-0.45, np.log10(10.**(-0.45)/1e3)+0.04, "1 Gyr", color="0.5", 
          ha="left", va="bottom", 
         rotation="25")
plt.plot(np.log10(xvec), np.log10(xvec/1e4), color="0.5", linestyle="--",
         linewidth=2)
plt.text(-0.45, np.log10(10.**(-0.45)/1e4)+0.04, "10 Gyr", color="0.5", 
          ha="left", va="bottom", 
          rotation="25")
plt.plot(np.log10(xvec), np.log10(xvec/1e5), color="0.5", linestyle="--",
         linewidth=2)
plt.text(-0.45, np.log10(10.**(-0.45)/1e5)+0.04, "100 Gyr", color="0.5", 
          ha="left", va="bottom", 
          rotation="25")
plt.plot(np.log10(xvec), np.log10(xvec/1e6), color="0.5", linestyle="--",
         linewidth=2)
plt.text(0.1, np.log10(10.**0.1/1e6)+0.04, "1000 Gyr", color="0.5", 
         ha="left", va="bottom", 
         rotation="25")

# Limits and labels
xlim([-0.5,2.5])
ylim([-6,0])
xlabel(r'$\log\,\Sigma_{\rm H_2}$ $[M_\odot\,{\rm pc}^{-2}]$')
ylabel(r'$\log\,\Sigma_{\rm SFR}$ $[M_\odot\,{\rm pc}^{-2}\,{\rm Myr}^{-1}]$')

# Legend
pbl = Rectangle((0,0), 1, 1, fc='b')
legend([pbl, p5, pgr, p3, p4, p1, p2], 
       ["Inner disks",
        "Inner disks (median)",
        "SMC (12 pc)",
        "SMC (200 pc)",
        "SMC (1 kpc)",
        "Rings (strong detection)",
        "Rings (weak detection)"],
       loc=4, numpoints=1)

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


#savefig('sflawh2plot.pdf')
