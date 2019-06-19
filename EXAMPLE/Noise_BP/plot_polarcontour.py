# -*- coding: utf-8 -*-
"""
Plot polar contour of beam power
18 June 2019
"""

#def BeampowerPolarplot(finame):

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import numpy as np
import sys
import math
import time
import h5py

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

#------Model Values-----#
finame = './dataset/BP_polarplot.h5'
colormap = cm.jet
A0amp = 0.1 # reference coefficient of magnitude of beam power to normalize
dBmin = -6 # minimum dB for plot
dBmax = 3 # maximum dB for plot
levelnum = 201; # level of color bar
whitecolor = (0.75, 0.75, 0.75)
freqvec = [0, 3, 6]
tstep = 2 # step of time for plot
#-----------------------#

print('#---------------------------#')
print('Plot Beam Power on Polar axis')
print('#---------------------------#')

f=h5py.File(finame,"r")

DLtimestamplist = f['info/DLtimestamplist'].value

#for i in range(len(DLtimestamplist)):
for i in [1, 2]:
    fvec = f['BeamPower/%s/fvec'%DLtimestamplist[i]].value
    fvecstamplist = f['BeamPower/%s/fvecstamplist'%DLtimestamplist[i]].value
    for fstamp in fvecstamplist[freqvec]:
    #for fstamp in ['f:0.20Hz']:
        cvec = f['BeamPower/%s/%s/cvec'%(DLtimestamplist[i], fstamp)].value * 1e-3 #[km/s]
        thetavec = np.deg2rad(f['BeamPower/%s/%s/thetavec'%(DLtimestamplist[i], fstamp)].value + 180) # to represent the location of source
        freq = f['BeamPower/%s/%s/f'%(DLtimestamplist[i], fstamp)].value
        time = f['BeamPower/%s/%s/time'%(DLtimestamplist[i], fstamp)].value
        
        # evaluate A0
        A0max=[]
        for tstamp in time[::tstep]:
            bp = np.transpose(f['BeamPower/%s/%s/%s/bp'%(DLtimestamplist[i], fstamp, tstamp)].value)
            A0max.append(np.amax(bp))

        A0 = A0amp * max(A0max)
        tcount = 0;
        for tstamp in time[::tstep]:
        #for tstamp in time[0:3]:
            tcount = tcount + 1
            bp = np.transpose(f['BeamPower/%s/%s/%s/bp'%(DLtimestamplist[i], fstamp, tstamp)].value)
            #bp to be dB
            bp_dB = 10*np.log10(bp/A0);
            #print(np.amax(bp_dB))
            r, theta = np.meshgrid(cvec, thetavec)
            #-- Plot... ------------------------------------------------
            with plt.rc_context({'ytick.color':whitecolor}):
                fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8, 8))
                levels = np.linspace(dBmin, dBmax, levelnum)
                cs = ax.contourf(theta, r, bp_dB, cmap=colormap, levels=levels, extend='both')
                # plot decoration
                ax.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])
                ax.set_yticks([2,3,4,5])
                ax.annotate('[km/s]', (1.00, 0.701), xycoords='axes fraction')
                ax.set_title('f: %4.2f [Hz] Time: %s'%(freq, tstamp), y=1.08)
                cbar = fig.colorbar(cs,  ticks=[-6, -4, -2, 0, 2, 4], orientation="horizontal", fraction=0.05, pad=0.1)
                cbar.set_label('Beam power [dB]')
                #plt.show()
                foname='./fig/beampower_%s_%s_%04d.png'%(DLtimestamplist[i], fstamp, tcount)
                plt.savefig(foname, dpi=100, facecolor='w', edgecolor='w')

