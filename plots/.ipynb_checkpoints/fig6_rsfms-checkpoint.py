import matplotlib
matplotlib.use('agg')

matplotlib.rcParams['font.size']=50

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
import h5py
from matplotlib.colors import LogNorm

h = 0.6774
fig_name = 'hist2d_sfr.pdf'

###################################################################
# data_dir = '../fixed/'

data_dir = '/orange/paul.torrey/jqi/TNG_resolved_scaling_relation/fixed_x/'
save_dir = '../'

# Load in Spaxel Data
file = h5py.File(data_dir + 'test_data2.hdf5', 'r')

file_list = sorted(glob.glob(data_dir + 'res/test_res_bin???.txt'))
names = ['mass090', 'mass095', 'mass100', 'mass105', 'mass110']
print file_list
print names
titles = ['$10^{9.0}\,<\,\mathrm{M}_{\star}\,<\,10^{9.5}\,\mathrm{M}_{\odot}$', \
        '$10^{9.5}\,<\,\mathrm{M}_{\star}\,<\,10^{10.0}\,\mathrm{M}_{\odot}$', \
        '$10^{10.0}\,<\,\mathrm{M}_{\star}\,<\,10^{10.5}\,\mathrm{M}_{\odot}$', \
        '$10^{10.5}\,<\,\mathrm{M}_{\star}\,<\,10^{11.0}\,\mathrm{M}_{\odot}$', \
        '$10^{11.0}\,<\,\mathrm{M}_{\star}\,<\,10^{15.0}\,\mathrm{M}_{\odot}$']

fig, axs = plt.subplots(2,2, figsize=(30.0,30.0), sharex=True, sharey=True)
axs = axs.flatten()
fig_name = 'hist2d_sfr_bins.pdf'

for k in xrange(len(names)-1):
    grp = file[names[k]]
    spx_mass  = np.array(grp['spx_stellar_mass'])
    spx_sfr   = np.array(grp['spx_sfr']         )
    spx_metal = np.array(grp['spx_metallicity'] )
    spx_gm    = np.array(grp['spx_gas_mass']    )
    
    spx_mass  = np.log10(spx_mass             ) + 6.0
    spx_sfr   = np.log10(spx_sfr * h          )
    spx_metal = np.log10(spx_metal * 0.35 / 16) + 12.0
    spx_gm    = np.log10(spx_gm               ) + 6.0

    mass, mass1, mass2, sfr, sfr1, sfr2, z, z1, z2 = np.loadtxt( file_list[k], unpack=True, usecols=(0,1,2,3,4,5,6,7,8))
    # mass, mass1, mass2, sfr, sfr1, sfr2, z, z1, z2 = np.loadtxt(
    #     data_dir + file_list[k], unpack=True, usecols=(0,1,2,3,4,5,6,7,8)
    # )

    mass  = np.log10(mass)  + 6.0
    mass1 = np.log10(mass1) + 6.0
    mass2 = np.log10(mass2) + 6.0
    
    sfr   = np.log10(sfr * h)
    sfr1  = np.log10(sfr1 * h)
    sfr2  = np.log10(sfr2 * h)

    z     = np.log10(z)
    z1    = np.log10(z1)
    z2    = np.log10(z2)

    sfr_l = -4.0
    sfr_h = 0.0
    sm_l  = 7.0
    sm_h  = 10.0

    nbins = 20
    sfr_sp = (sfr_h - sfr_l) / nbins
    sm_sp = (sm_h - sm_l) / nbins

    m2 = np.zeros((nbins, nbins))
    m1 = [[[] for _ in range(nbins) ] for _ in range(nbins)]

    for i in xrange(len(spx_mass)):
        if spx_sfr[i] < sfr_l: continue
        y = nbins - 1 - int( (spx_sfr[i] - sfr_l) / sfr_sp )
        x = int( (spx_mass[i] - sm_l) / sm_sp )
        if (x >= 0) and (x < nbins) and (y >= 0) and (y < nbins):
                m1[y][x].append(spx_metal[i])
                m2[y][x] += 1
            
    m3 = np.ones((nbins, nbins))
    for i in range(nbins):
        for j in range(nbins):
            if m2[i][j] > 10:
                m3[i][j] = m2[i][j]
            else:
                m3[i][j] = None
                
    # Plot Figure
	ax = axs[k]

    ax.tick_params(labelsize=50)
	# mpb = plt.hist2d(spx_mass, spx_sfr, bins=20, norm=LogNorm(), cmap='coolwarm', cmin=10, cmax=5000)
    mpb = ax.imshow(m3, extent=(sm_l, sm_h, sfr_l, sfr_h), norm=LogNorm(), interpolation='none', vmin=10, vmax=3000, cmap='coolwarm', aspect='auto')
    
    ax.plot(mass, sfr, 'black', lw=5, label='IllustrisTNG')
    ax.plot(mass, sfr1, 'black', ls='--', lw=3)
    ax.plot(mass, sfr2, 'black', ls='--', lw=3)
    
    if (k == 0 or k==2):
        ax.set_ylabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2} \, \mathrm{yr}^{-1})$', fontsize=50)
    if (k == 2 or k==3):
        ax.set_xlabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{\star}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=50)

	ax.tick_params(labelsize=50)
    ax.set_xlim([7.01, 9.99])
    ax.set_ylim([-4.0, -0.01])
    
    ax.text( 0.45,0.9, titles[k] ,transform=ax.transAxes, fontsize=50 )
    
ax_cbar = fig.add_axes([0.2, 1.02, 0.6, 0.02])

cb = plt.colorbar(mpb, cax=ax_cbar,orientation='horizontal',aspect=10,shrink=0.5)
# cb.set_label(r'${\rm Radius~(kpc)}$',fontsize=50,labelpad=100)

cb.ax.text( 0.5,1.27, r'${\rm Count}$' , transform=cb.ax.transAxes, fontsize=50, ha='center' )

# cb.ax.xaxis.set_ticks_position('top')
cb.ax.xaxis.set_label_position('top')

cb.ax.tick_params(labelsize=50)

plt.tight_layout()

plt.subplots_adjust(wspace=0.01, hspace=0.01)
plt.savefig(save_dir+fig_name,bbox_inches='tight')

file.close()	
 #    fig = plt.figure(figsize=(3.0, 15.0))
 #    ax = fig.add_axes([0.10, 0.10, 0.82, 0.82])
 #    ax.remove()
 #    cbaxes = fig.add_axes([0.02, 0.10, 0.22, 0.82])
 #    cbar = plt.colorbar(mpb, cax=cbaxes, pad=0.03, fraction=0.08, ticks=[10**1, 10**2, 10**3])
 #    cbar.set_label('Count', fontsize=50, labelpad=1.1)
 #    cbar.ax.tick_params(labelsize=50)
 #    cbar.ax.xaxis.set_ticks_position('top')
 #    cbar.ax.xaxis.set_label_position('top')
	# plt.savefig('cnt_cb.png')

file.close()
