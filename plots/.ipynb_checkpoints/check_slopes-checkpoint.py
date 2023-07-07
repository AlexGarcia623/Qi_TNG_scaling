import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt

import h5py
import numpy as np
from scipy.interpolate import UnivariateSpline

data_dir = '/orange/paul.torrey/jqi/TNG_resolved_scaling_relation/fit_all/'
save_dir = '../'


file = h5py.File(data_dir +  'values.hdf5', 'r')

print file.keys()
gal_bin = [9.0, 9.5, 10.0, 10.5, 11.0, 15.0]
intercepts = [[], [], [], [], []]
slopes     = [[], [], [], [], []]
scatter    = [[], [], [], [], []]

norms      = [[], [], [], [], []]

itp = file['intercept']
s   = file['slopes']
sca = file['scatter']

print sca.keys()

for i in xrange(len(slopes)):
	intercepts[i] = np.array(itp[itp.keys()[i]])
	slopes[i]     = np.array(s[s.keys()[i]])
	scatter[i]    = np.array(sca[sca.keys()[i]])

all_slp = np.array([])
all_sca = np.array([])
all_nom = np.array([])

sm_standard = 8.0

for i in xrange(len(slopes)):
        for j in xrange(len(slopes[i])):
                this_norm = sm_standard * slopes[i][j] + intercepts[i][j]
		norms[i].append(this_norm)

for i in xrange(len(slopes) - 1):
	all_slp = np.concatenate((all_slp, np.array(slopes[i])))
	all_sca = np.concatenate((all_sca, np.array(scatter[i])))
	all_nom = np.concatenate((all_nom, np.array(norms[i])))

colors = ['red', 'orange', 'green', 'blue']

print len(all_slp), len(all_sca), len(all_nom)
for i in xrange(len(slopes)):
	print len(intercepts[i]), len(slopes[i]), len(scatter[i])

ps = [0, 0, 0, 0, 0]
xs = [0, 0, 0, 0, 0]
p, x, patches = plt.hist(all_nom, 30, alpha=0.75)
for i in xrange(len(norms)):
        ps[i], xs[i], patches = plt.hist(np.array(norms[i]), 30, alpha=0.75)

fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.13, 0.10, 0.80, 0.80])
avg_nom = np.mean(all_nom)
med_nom = np.median(all_nom)
x = x[:-1] + (x[1] - x[0]) / 2.0
f = UnivariateSpline(x, p, s=30)
ax.plot(x, f(x) / float(len(all_nom)), lw=7, color='black')

for i in xrange(len(slopes) - 1):
        x = xs[i][:-1] + (xs[i][1] - xs[i][0]) / 2.0
        f = UnivariateSpline(x, ps[i], s=30)
        ax.plot(x, f(x) / float(len(slopes[i])), label= str(gal_bin[i]) + ' - ' + str(gal_bin[i+1]), lw=5, color=colors[i])

ax.set_xlabel('Normalization of $\mathrm{\Sigma}_{\mathrm{SFR}}$ at $\mathrm{\Sigma}_{\star}\,=\,10^{8.0}\,\mathrm{M}_{\odot}\,\mathrm{kpc}^{-2}$', fontsize=50)
ax.set_ylabel('Probability', fontsize=50)
ax.tick_params(labelsize=50)
ax.set_title("Mean={:2.4f}, Median={:2.4f}".format(avg_nom, med_nom), fontsize=50)
#plt.legend(fancybox=True, framealpha=0.5, fontsize=40)
ax.set_xlim([-4.0, -1.0])
ax.set_ylim([0.01, 0.3])
plt.savefig(save_dir + 'all_norms.pdf')

ps = [0, 0, 0, 0, 0]
xs = [0, 0, 0, 0, 0]
p, x, patches = plt.hist(all_slp, 30, alpha=0.75)
for i in xrange(len(slopes)):
	ps[i], xs[i], patches = plt.hist(np.array(slopes[i]), 30, alpha=0.75)

fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.13, 0.10, 0.80, 0.80])
avg_slp = np.mean(all_slp)
med_slp = np.median(all_slp)
x = x[:-1] + (x[1] - x[0]) / 2.0
f = UnivariateSpline(x, p, s=30)
ax.plot(x, f(x) / float(len(all_slp)), lw=7, color='black')

labels = [
	r'$10^{9.0} <M_*<10^{9.5} M_\odot$',
	r'$10^{9.5} <M_*<10^{10.0}M_\odot$',
	r'$10^{10.0}<M_*<10^{10.5}M_\odot$',
	r'$10^{10.5}<M_*<10^{11.0}M_\odot$',
	'all'
]

for i in xrange(len(slopes) - 1):
	x = xs[i][:-1] + (xs[i][1] - xs[i][0]) / 2.0
	f = UnivariateSpline(x, ps[i], s=30)

	ax.plot(x, f(x) / float(len(slopes[i])), label= str(gal_bin[i]) + ' - ' + str(gal_bin[i+1]), lw=5, color=colors[i])

ax.set_xlabel('Slope', fontsize=50)
ax.set_ylabel('Probability', fontsize=50)
ax.tick_params(labelsize=50)
ax.set_title("Mean={:2.4f}, Median={:2.4f}".format(avg_slp, med_slp), fontsize=50)
#plt.legend(fancybox=True, framealpha=0.5, fontsize=40)
ax.set_xlim([-1.5, 3.0])
ax.set_ylim([0.01, 0.3])
plt.savefig(save_dir + 'all_slopes.pdf')



p, x, patches = plt.hist(all_sca, 30, alpha=0.75)
ps = [0, 0, 0, 0, 0]
xs = [0, 0, 0, 0, 0]
for i in xrange(len(scatter)):
        ps[i], xs[i], patches = plt.hist(np.array(scatter[i]), 30, alpha=0.75)

fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.13, 0.10, 0.80, 0.80])
avg_sca = np.mean(all_sca)
med_sca = np.median(all_sca)
x = x[:-1] + (x[1] - x[0]) / 2.0
f = UnivariateSpline(x, p, s=30)
ax.plot(x, f(x) / float(len(all_sca)), lw=7, color='black', label='All Galaxies')

for i in xrange(len(slopes) - 1):
        x = xs[i][:-1] + (xs[i][1] - xs[i][0]) / 2.0
        f = UnivariateSpline(x, ps[i], s=30)

	ax.plot(x, f(x) / float(len(slopes[i])), label=labels[i], lw=5, color=colors[i])

#        ax.plot(x, f(x) / float(len(slopes[i])), label= str(gal_bin[i]) + ' - ' + str(gal_bin[i+1]), lw=5, color=colors[i])

ax.set_xlabel('Scatter / dex', fontsize=50)
ax.set_ylabel('Probability', fontsize=50)
ax.tick_params(labelsize=50)
ax.set_title("Mean={:2.4f}, Median={:2.4f}".format(avg_sca, med_sca), fontsize=50)
plt.legend(fancybox=True, framealpha=0.5, fontsize=50)
ax.set_xlim([0.0, 0.8])
ax.set_ylim([0.01, 0.3])
plt.savefig(save_dir + 'all_scatters.pdf')
