import astropy.table as tab
from astropy.io import ascii
import glob
import numpy as np

path  = '/home/vikram/Dropbox/power_spectrum_idl_cleaned_up/danforth_data_MW_format_new/modified_MW_format'
file = glob.glob(path + '/*.dat')

# making 40 angstrom chuncks
chuncks = np.arange(17)*40+1120
masked_wave = []
for c in chuncks:
    print('on chunck', c, 'to', c+40)
    wave = []
    for f in file:
        data = ascii.read(f, data_start=1)
        wave_new = data['col1'][data['col5']==7]
        wave_new  = wave_new[(wave_new>=c) & (wave_new<(c+40))]
        if len(wave_new)> len(wave):
            wave = wave_new
    print('size', len(wave))
    if len(wave)>0:
        masked_wave = masked_wave + list(wave)

wave_table = tab.Table([masked_wave], names=['mask'])
wave_table.write('milky_way_masks.fits', overwrite = True)