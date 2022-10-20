import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import astropy.table as tab
import numpy as np
import scipy.interpolate as interp
from scipy.signal import savgol_filter
# example
#yhat = savgol_filter(y, 51, 3) # window size 51, polynomial order 3


def readfile(filename):
   ps0=tab.Table.read(filename, format='ascii.no_header',names=['z','k','kpk','errkpk','N'],guess=False)
 #  ps0=ps0[ps0['k']>9e-4]
 #  ps0=ps0[ps0['k']<0.12]
   ps0=ps0[ps0['N']>=4]
   x=ps0['k']
   y=ps0['kpk']
   err=ps0['errkpk']
   z=ps0['z'][0]
   #label='z = {:.2f} ({:.3f} - {:.2f})'.format(z, binz[0], binz[1])
   return x, y, err, z

# plotting
fig_name = 'power_hsla_z01_scale_cont.pdf'
font = {'family': 'serif', 'weight': 'normal', 'size': 11}
plt.rc('font', **font)
fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, gridspec_kw = {'height_ratios':[3, 1]},figsize=(6, 7) )

# Read the power-spectrum file
binz = [0.06, 0.16]
filename='/home/vikram/power_data/z0Power_metals_masked_z{:.2f}-{:.2f}_wave1050-1180.txt'.format(binz[0], binz[1])
# same folder is there in igm server
x, y, err, z=readfile(filename)
label=r'$\bar {\rm z}$'+'={:.2f} (Khaire et al. 2019)'.format(z)

ax1.errorbar(x, y, err, label=label, marker='.', markersize=12, ls='', color='green', capsize=5, elinewidth=2, zorder = 10, alpha =0.75)

filename='/home/vikram/output_power/hsla_power/high_SN/Power_metals_masked_z{:.2f}-{:.2f}_wave1050-1180.txt'.format(binz[0], binz[1])
# same folder is there in igm server
x, y, err, z=readfile(filename)
label=r'$\bar {\rm z}$'+'={:.2f} HSLA high SN'.format(z)

ax1.errorbar(x, y, err, label=label, marker='.', markersize=12, ls='', color='red', capsize=5, elinewidth=2, zorder = 10)

"""
filename='/home/vikram/output_power/scale_cont/Power_metals_masked_z{:.2f}-{:.2f}_wave1050-1180.txt'.format(binz[0], binz[1])
# same folder is there in igm server
x, y, err, z=readfile(filename)
label=r'$\bar {\rm z}$'+'={:.2f} cont*1.05'.format(z)

ax1.errorbar(x, y, err, label=label, marker='.', markersize=12, ls='', color='k', capsize=5, elinewidth=2, zorder = 10, alpha =0.5)

"""


filename='/home/vikram/output_power/hsla_power/medium_3.1_quality/Power_metals_masked_z{:.2f}-{:.2f}_wave1050-1180.txt'.format(binz[0], binz[1])
# same folder is there in igm server
x, y, err, z=readfile(filename)
label=r'$\bar {\rm z}$'+'={:.2f} HSLA Medium SN'.format(z)
ax1.errorbar(x, y, err, label=label, marker='.', markersize=12, ls='', color='orange', capsize=5, elinewidth=2, zorder = 10, alpha=0.4)

ax2.set_ylim (-12, 12)
ax2.set_ylabel ('% difference')
#ax2.scatter(k3, r3, color='g', label='', s=20)
#ax2.hlines(y=0.0, xmin=1e-10, xmax=1e10)
#ax2.plot([1e-10, 1e10], [0.0, 0.0], linestyle='--',dashes=(5, 6), color='k')
ax2.hlines(y=0, xmin=1e-10, xmax=1e10, linestyles = '--', linewidth=1, color= 'blue')
ax2.hlines(y=5, xmin=1e-10, xmax=1e10, linestyles=':', color='b')
ax2.hlines(y=-5, xmin=1e-10, xmax=1e10, linestyles=':', color='b')

ax1.set_yscale('log', nonposy='clip')
ax1.set_ylim(7e-5, 0.015)
ax1.legend(loc='best', fontsize = 11)
ax1.set_ylabel(r'k P(k) ${\rm \, /\, \pi}$')
ax2.set_xlabel(r'k [s/km]')
ax2.tick_params(axis='x', which='major', pad=10)  # to keep a gap between labels and line

for ax in (ax1, ax2):
   # decorating the plot
   ax.set_xlim(4e-4, 0.25)
   ax.set_xscale('log')
   ax.tick_params(direction='in', length=5, width=1.5)
   ax.tick_params(direction='in', which='minor', length=3.5, width=1.5)
   ax.xaxis.set_ticks_position('both')
   ax.yaxis.set_ticks_position('both')
   for axis in ['top', 'bottom', 'left', 'right']:
      ax.spines[axis].set_linewidth(1.5)


fig.tight_layout(h_pad=1)
fig.savefig(fig_name, bbox_inches='tight')
