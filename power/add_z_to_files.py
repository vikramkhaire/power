import astropy.table as tab
import numpy as np


masterfile = '/home/vikram/Dropbox/power_spectrum_idl_cleaned_up/danforth_data_MW_format_new/masterfile.txt'

forward_file ='/home/vikram/output_power/test/danforth_forward_models_min_SN_10_Gamma_0.09000_Nran_010000_seed_42.fits'

data = tab.Table.read(forward_file)
obj_fwd = data['corresponding_data_obj']

obj_obs, redshift, sn, flag, res, lpin = zip(
    *np.genfromtxt(masterfile, comments='#', dtype=['<U30', np.float, np.float, np.int, np.float, 'U3']))

z_fwd = []

for i in range(len(data)):
    ind = obj_obs.index(obj_fwd[i])
    z_fwd.append(redshift[ind])

data.add_column(z_fwd, name = 'zqso')
print('adding zqso column')
data.write(forward_file, overwrite = True)
