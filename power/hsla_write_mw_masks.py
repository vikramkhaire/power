import astropy.table as tab
import numpy as np
import glob


data_path = '/home/vikram/output_power/hsla'

files =  glob.glob(data_path + '/*final*.fits')

mw_masks_table = tab.Table.read('milky_way_masks.fits')
mw_wave = mw_masks_table['mask']

for f in files:
    data = tab.Table.read(f)
    wave = data['WAVE']
    delta_wave = wave[1]-wave[0]
    new_masks  =[]
    for l in wave:
        diff = np.abs(mw_wave-l)
        min_index = np.argmin(diff)
        if diff[min_index]<2*delta_wave:
            new_masks.append(7)
            #print(l, mw_wave[min_index], '7')
        else:
            new_masks.append(0)
            #print(l, mw_wave[min_index], '0')
    data.add_column(new_masks, name = 'mask')
    data.write(f, overwrite= True)
    print('written', f)

