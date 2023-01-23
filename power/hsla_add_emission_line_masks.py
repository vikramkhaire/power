import astropy.table as tab
import numpy as np
import glob
from astropy.io import ascii

# assumes that the milky way metal masks are added (i.e a column name mask exists)

path  = '/home/vikram/Dropbox/power_spectrum_idl_cleaned_up/danforth_data_MW_format_new/modified_MW_format'
file = glob.glob(path + '/*.dat')

"""
low = []
high = []
for f in file:
    data = ascii.read(f, data_start=1)
    data_red =data[data['col1']> 1230]
    wave_red= data_red['col1'][data_red['col5']==6]
    if len(wave_red)>2:
        print(min(wave_red), max(wave_red))
        low.append(min(wave_red))
        high.append(max(wave_red))

print(min(low), ': minimum emission line mask wavelength')
print(max(high), ': max --')
"""

# results of 1300.985107421875 : minimum emission line mask wavelength
#1300.985107421875 : minimum emission line mask wavelength
#1307.4949951171875 : max --


data_path = '/home/vikram/output_power/data'

files =  glob.glob(data_path + '/*final*.fits')


for f in files:
    data = tab.Table.read(f)
    wave = data['WAVE']
    data['mask'][(data['WAVE']>1300.98)& (data['WAVE']<1307.5)]=6
    print(len(data['mask'][data['mask']==6]), ' :size', f)

    data.write(f, overwrite= True)
    print('written', f)


import astropy.table as tab
def replace_mask(spectrum_file, maskvalue, wavelength_range= [1220, 1330]):
    """
    :param spectrum_file: fits file in which you are replacing mask
    :param maskvalue: value of the mask e.g 4
    :param wavelength_range: lower and higer wavelength in a list of two elements where you want to add maskvalue
    :return:
    """
    data  = tab.Table.read(spectrum_file)
    data['mask'][(data['WAVE']>wavelength_range[0])& (data['WAVE']<wavelength_range[1])]=maskvalue
    data.write(f, overwrite= True)
    print('written', f)

    return