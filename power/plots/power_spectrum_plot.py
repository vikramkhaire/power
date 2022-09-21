import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import astropy.table as tab
import numpy as np
import scipy.interpolate as interp
from scipy.signal import savgol_filter
# example
#yhat = savgol_filter(y, 51, 3) # window size 51, polynomial order 3

path = '/home/vikram/output_power/test'
LS_power_file = path + '/' + 'power_forward_tng_Gamma_0.09000_Nran_010000_seed_42.txt'
normal_power_file  =  path + '/' + 'new_power_igm_IllustrisTNG_z01_rudi_G0.09.fits'

d = tab.Table. read(LS_power_file, format = 'ascii')
p = tab.Table.read(normal_power_file)


# Read the power-spectrum file
binz = [0.06, 0.16]
filename='/home/vikram/power_data/z0Power_metals_masked_z{:.2f}-{:.2f}_wave1050-1180.txt'.format(binz[0], binz[1])
# same folder is there in igm server
#x, y, err, z=readfile(filename)
#label=r'$\bar {\rm z}$'+'={:.2f} (Khaire et al. 2019)'.format(z)


# plotting
fig_name = 'power_test.pdf'


font = {'family': 'serif', 'weight': 'normal', 'size': 11}
plt.rc('font', **font)
fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, gridspec_kw = {'height_ratios':[3, 1]},figsize=(6, 7) )


label = 'Lomb-scargle'
ax1.plot(d['LS K-modes'], d['LS Power(noise and wind corr)'],  label = label, color = 'b', linewidth=2.5)
label = 'normal'
ax1.plot(p['k_v'], p['pkp_v'], label = label, color = 'r', linestyle = '--', linewidth=2.5, dashes=(4,3))




logkpki=interp.interp1d(np.log10(p['k_v']), np.log10(p['pkp_v']),fill_value='extrapolate')
new_pkp_normal = 10**(logkpki(np.log10(d['LS K-modes'][~np.isnan(d['LS K-modes'])]))) # at x values of old ill
r_ill_by_tng=(d['LS Power(noise and wind corr)'][~np.isnan(d['LS K-modes'])] - new_pkp_normal)/ new_pkp_normal

ax2.plot(d['LS K-modes'][~np.isnan(d['LS K-modes'])], r_ill_by_tng*100, color = 'k')

ax2.set_ylim (-12, 12)
ax2.set_ylabel ('% difference')
#ax2.scatter(k3, r3, color='g', label='', s=20)
#ax2.hlines(y=0.0, xmin=1e-10, xmax=1e10)
#ax2.plot([1e-10, 1e10], [0.0, 0.0], linestyle='--',dashes=(5, 6), color='k')
ax2.hlines(y=0, xmin=1e-10, xmax=1e10, linestyles = '--', linewidth=1, color= 'blue')
ax2.hlines(y=5, xmin=1e-10, xmax=1e10, linestyles=':', color='b')
ax2.hlines(y=-5, xmin=1e-10, xmax=1e10, linestyles=':', color='b')

ax1.set_yscale('log', nonposy='clip')
ax1.set_ylim(7e-5, 0.01)
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
