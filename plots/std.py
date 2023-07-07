import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
import h5py
from matplotlib.colors import LogNorm
h = 0.6774

###################################################################
# data_dir = '../fixed/'
data_dir = '/orange/paul.torrey/jqi/Share/AlexG/TNG_resolved_scaling_relation/L35n2160TNG/'
save_dir = '../'

ASPECT = 0.75

gm, gm1, gm2, sfr, sfr1, sfr2 = \
        np.loadtxt(data_dir + 'res/res_ks110.txt', unpack=True, usecols=(0,1,2,3,4,5))

file = h5py.File(data_dir + 'spx_data.hdf5', 'r')

spx_gm  = np.log10(np.array(file['ks/spx_gas_mass'])) + 6.0
spx_sfr = np.log10(np.array(file['ks/spx_sfr']) * h)
spx_sm  = np.log10(np.array(file['ks/spx_stellar_mass'])) + 6.0
spx_fg  = np.log10(np.array(file['ks/spx_gas_ratio']))
spx_z   = np.log10(np.array(file['ks/spx_gas_metallicity']))
spx_rad =          np.array(file['tot/spx_radius'])
print len(spx_rad), len(spx_gm)
file.close()

#################################################################
sfr_l = -4.0
sfr_h = 0.0
gm_l  = 6.5
gm_h  = 9.0

nbins  = 20
sfr_sp = (sfr_h - sfr_l) / nbins
gm_sp  = (gm_h - gm_l) / nbins

#m1 = np.zeros((nbins, nbins))
m2 = np.zeros((nbins, nbins))

m1 = [[[] for _ in range(nbins) ] for _ in range(nbins)]

for i in xrange(len(spx_gm)):
	y = nbins - 1 - int( (spx_sfr[i] - sfr_l) / sfr_sp )
	x = int( (spx_gm[i] - gm_l) / gm_sp )
	if (x >= 0) and (x < nbins) and (y >= 0) and (y < nbins):
		#m1[y][x] += spx_sfr[i]
		m1[y][x].append(spx_sm[i])
		m2[y][x] += 1

#m3 = np.divide(m1, m2, out=np.zeros_like(m2), where=m2>5)
m3 = np.zeros((nbins, nbins))
for i in range(nbins):
    for j in range(nbins):
        if m2[i][j] > 10:
            # m3[i][j] = np.median(m1[i][j])
            m3[i][j] = m2[i][j]
        else:
            m3[i][j] = None


fig, axs = plt.subplots(1,2,figsize=(25.0, 15.0),gridspec_kw={'height_ratios': [1],'width_ratios':[1,1]})

axs = axs.flatten()

############################################################### START PANEL 1 ###############################################################
ax = axs[0]
histogram = ax.imshow(m3, extent=(gm_l, gm_h, sfr_l, sfr_h), interpolation='none', cmap='coolwarm', norm=LogNorm(), aspect=ASPECT)
# ax.hist2d(spx_gm, spx_sfr, bins=20, norm=LogNorm(), cmap='coolwarm',vmin=10)
gm  = np.log10(gm) + 6.0
sfr = np.log10(sfr * h)
ax.plot(gm, sfr, 'black', lw=5, label='All Galaxies')

file_list = sorted(glob.glob(data_dir + '/res/ks1_bin*'))

mass_bin = [9.0, 9.5, 10.0, 10.5, 11.0]
colors = ['red', 'orange', 'green', 'blue']

for i in xrange(len(file_list) - 1):
        x, y = np.loadtxt(file_list[i], unpack=True, usecols=(0,3))
        x = np.log10(x) + 6.0
        y = np.log10(y * h)
	if (i < len(file_list) - 1):
	        ax.plot(x, y, color=colors[i], label= str(mass_bin[i]) + ' - ' + str(mass_bin[i+1]) + r' $\log M_\odot$' , lw=3)

ax.legend(fancybox=True, framealpha=0.5, fontsize=35)

ax.set_xlabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{gas}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=50)
ax.set_ylabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2} \, \mathrm{yr}^{-1})$', fontsize=50)
ax.tick_params(labelsize=50)
ax.set_xlim([6.01, 9.0])
ax.set_ylim([-4.0, -0.01])

# plt.savefig(save_dir+'ks_count.pdf')

############################################################### START PANEL 2 ###############################################################
ax = axs[1]
sm_l = 7.0
sm_h = 10.0
fg_l = -3.0
fg_h = 0.0

nbins = 20
sm_sp = (sm_h - sm_l) / nbins
fg_sp = (fg_h - fg_l) / nbins

#m1 = np.zeros((nbins, nbins))
m2 = np.zeros((nbins, nbins))

m1 = [[[] for _ in range(nbins) ] for _ in range(nbins)]

for i in xrange(len(spx_sm)):
        y = nbins - 1 - int( (spx_fg[i] - fg_l) / fg_sp )
        x = int( (spx_sm[i] - sm_l) / sm_sp )
        if (x >= 0) and (x < nbins) and (y >= 0) and (y < nbins):
                #m1[y][x] += spx_sfr[i]
                m1[y][x].append(spx_sfr[i])
                m2[y][x] += 1

#m3 = np.divide(m1, m2, out=np.zeros_like(m2), where=m2>5)
m3 = np.zeros((nbins, nbins))
for i in range(nbins):
    for j in range(nbins):
        if m2[i][j] > 10:
            # m3[i][j] = np.median(m1[i][j])
            m3[i][j] = m2[i][j]
        else: 
            m3[i][j] = None
            

ax.imshow(m3, extent=(sm_l, sm_h, fg_l, fg_h), interpolation='none', cmap='coolwarm', norm=LogNorm(), aspect=ASPECT*1.33)
# histogram = ax.hist2d(spx_sm, spx_fg, bins=20, norm=LogNorm(), cmap='coolwarm', vmin=10)

sm, sm1, sm2, fg, fg1, fg2 = \
        np.loadtxt(data_dir + 'res/res_ks210.txt', unpack=True, usecols=(0,1,2,3,4,5))

sm = np.log10(sm) + 6.0
fg = np.log10(fg)

ax.plot(sm, fg, 'black', lw=5)

# cbar = plt.colorbar(pad=0.03, orientation='horizontal', fraction=0.08)
# cbar.set_label('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2}\, \mathrm{yr}^{-1})$', fontsize=50, rotation=90, labelpad=1.1)
# cbar.ax.tick_params(labelsize=40)

file_list = sorted(glob.glob(data_dir + 'res/ks2_bin*'))

mass_bin = [9.0, 9.5, 10.0, 10.5, 11.0]
colors = ['red', 'orange', 'green', 'blue']

for i in xrange(len(file_list) - 1):
        x, y = np.loadtxt(file_list[i], unpack=True, usecols=(0,3))
        x = np.log10(x) + 6.0
        y = np.log10(y)
        if (i < len(file_list) - 1):
                ax.plot(x, y, color=colors[i], label= str(mass_bin[i]) + ' - ' + str(mass_bin[i+1]), lw=3)

ax.set_xlabel('$\mathrm{log}(\mathrm{\Sigma}_{\star} / M_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=50)
ax.set_ylabel('$\mathrm{log}(\mathrm{f}_{\mathrm{gas}}) $', fontsize=50)
ax.tick_params(labelsize=50)
ax.set_xlim([7.01, 10.0])
ax.set_ylim([-3.0, 0.0])

cbaxes = fig.add_axes([0.175, 0.9, 0.7, 0.05])#[0.2, 1.02, 0.6, 0.02])
# cbaxes = fig.add_subplot(123)
cbar = plt.colorbar(histogram, cax=cbaxes, pad=0.03, orientation='horizontal', fraction=0.08)
# cbar.set_label('$\mathrm{Count}$', fontsize=50, labelpad=1.2)
cbar.ax.text( 0.5,1.2, r'${\rm Count}$', transform=cbar.ax.transAxes, fontsize=50, ha='center') 
cbar.ax.tick_params(labelsize=50)
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.ax.xaxis.set_label_position('top')

plt.tight_layout()

plt.savefig(save_dir+'ks_fgas_count.pdf')



#################################################################
sfr_l = -4.0
sfr_h = 0.0
gm_l  = 6.5
gm_h  = 9.0

nbins  = 20
sfr_sp = (sfr_h - sfr_l) / nbins
gm_sp  = (gm_h - gm_l) / nbins

#m1 = np.zeros((nbins, nbins))
m2 = np.zeros((nbins, nbins))

m1 = [[[] for _ in range(nbins) ] for _ in range(nbins)]

for i in xrange(len(spx_gm)):
	y = nbins - 1 - int( (spx_sfr[i] - sfr_l) / sfr_sp )
	x = int( (spx_gm[i] - gm_l) / gm_sp )
	if (x >= 0) and (x < nbins) and (y >= 0) and (y < nbins):
		#m1[y][x] += spx_sfr[i]
		m1[y][x].append(spx_rad[i])
		m2[y][x] += 1

#m3 = np.divide(m1, m2, out=np.zeros_like(m2), where=m2>5)
m3 = np.zeros((nbins, nbins))
for i in range(nbins):
    for j in range(nbins):
        if m2[i][j] > 10:
            m3[i][j] = np.median(m1[i][j])
            # m3[i][j] = m2[i][j]
        else:
            m3[i][j] = None


fig, axs = plt.subplots(1,2,figsize=(25.0, 15.0),gridspec_kw={'height_ratios': [1],'width_ratios':[1,1]})

axs = axs.flatten()

############################################################### START PANEL 3 ###############################################################
ax = axs[0]
histogram = ax.imshow(m3, extent=(gm_l, gm_h, sfr_l, sfr_h), interpolation='none', cmap='BrBG', aspect=ASPECT)
# ax.hist2d(spx_gm, spx_sfr, bins=20, norm=LogNorm(), cmap='coolwarm',vmin=10)
gm  = np.log10(gm) + 6.0
sfr = np.log10(sfr * h)
ax.plot(gm, sfr, 'black', lw=5, label='All Galaxies')

file_list = sorted(glob.glob(data_dir + '/res/ks1_bin*'))

mass_bin = [9.0, 9.5, 10.0, 10.5, 11.0]
colors = ['red', 'orange', 'green', 'blue']

for i in xrange(len(file_list) - 1):
        x, y = np.loadtxt(file_list[i], unpack=True, usecols=(0,3))
        x = np.log10(x) + 6.0
        y = np.log10(y * h)
	if (i < len(file_list) - 1):
	        ax.plot(x, y, color=colors[i], label= str(mass_bin[i]) + ' - ' + str(mass_bin[i+1]) + r' $\log M_\odot$' , lw=3)

ax.legend(fancybox=True, framealpha=0.5, fontsize=35)

ax.set_xlabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{gas}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=50)
ax.set_ylabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2} \, \mathrm{yr}^{-1})$', fontsize=50)
ax.tick_params(labelsize=50)
ax.set_xlim([6.01, 9.0])
ax.set_ylim([-4.0, -0.01])

# plt.savefig(save_dir+'ks_count.pdf')

############################################################### START PANEL 4 ###############################################################
ax = axs[1]
sm_l = 7.0
sm_h = 10.0
fg_l = -3.0
fg_h = 0.0

nbins = 20
sm_sp = (sm_h - sm_l) / nbins
fg_sp = (fg_h - fg_l) / nbins

#m1 = np.zeros((nbins, nbins))
m2 = np.zeros((nbins, nbins))

m1 = [[[] for _ in range(nbins) ] for _ in range(nbins)]

for i in xrange(len(spx_sm)):
        y = nbins - 1 - int( (spx_fg[i] - fg_l) / fg_sp )
        x = int( (spx_sm[i] - sm_l) / sm_sp )
        if (x >= 0) and (x < nbins) and (y >= 0) and (y < nbins):
                #m1[y][x] += spx_sfr[i]
                m1[y][x].append(spx_rad[i])
                m2[y][x] += 1

#m3 = np.divide(m1, m2, out=np.zeros_like(m2), where=m2>5)
m3 = np.zeros((nbins, nbins))
for i in range(nbins):
    for j in range(nbins):
        if m2[i][j] > 10:
            m3[i][j] = np.median(m1[i][j])
            # m3[i][j] = m2[i][j]
        else: 
            m3[i][j] = None
            

ax.imshow(m3, extent=(sm_l, sm_h, fg_l, fg_h), interpolation='none', cmap='BrBG', aspect=ASPECT*1.33)
# histogram = ax.hist2d(spx_sm, spx_fg, bins=20, norm=LogNorm(), cmap='coolwarm', vmin=10)

sm, sm1, sm2, fg, fg1, fg2 = \
        np.loadtxt(data_dir + 'res/res_ks210.txt', unpack=True, usecols=(0,1,2,3,4,5))

sm = np.log10(sm) + 6.0
fg = np.log10(fg)

ax.plot(sm, fg, 'black', lw=5)

# cbar = plt.colorbar(pad=0.03, orientation='horizontal', fraction=0.08)
# cbar.set_label('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2}\, \mathrm{yr}^{-1})$', fontsize=50, rotation=90, labelpad=1.1)
# cbar.ax.tick_params(labelsize=40)

file_list = sorted(glob.glob(data_dir + 'res/ks2_bin*'))

mass_bin = [9.0, 9.5, 10.0, 10.5, 11.0]
colors = ['red', 'orange', 'green', 'blue']

for i in xrange(len(file_list) - 1):
        x, y = np.loadtxt(file_list[i], unpack=True, usecols=(0,3))
        x = np.log10(x) + 6.0
        y = np.log10(y)
        if (i < len(file_list) - 1):
                ax.plot(x, y, color=colors[i], label= str(mass_bin[i]) + ' - ' + str(mass_bin[i+1]), lw=3)

ax.set_xlabel('$\mathrm{log}(\mathrm{\Sigma}_{\star} / M_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=50)
ax.set_ylabel('$\mathrm{log}(\mathrm{f}_{\mathrm{gas}}) $', fontsize=50)
ax.tick_params(labelsize=50)
ax.set_xlim([7.01, 10.0])
ax.set_ylim([-3.0, 0.0])

cbaxes = fig.add_axes([0.175, 0.9, 0.7, 0.05])#[0.2, 1.02, 0.6, 0.02])
# cbaxes = fig.add_subplot(123)
cbar = plt.colorbar(histogram, cax=cbaxes, pad=0.03, orientation='horizontal', fraction=0.08)
# cbar.set_label('$\mathrm{Count}$', fontsize=50, labelpad=1.2)
cbar.ax.text( 0.5,1.2, r'${\rm Radius/kpc}$', transform=cbar.ax.transAxes, fontsize=50, ha='center') 
cbar.ax.tick_params(labelsize=50)
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.ax.xaxis.set_label_position('top')

plt.tight_layout()

plt.savefig(save_dir+'ks_fgas_rad.pdf')