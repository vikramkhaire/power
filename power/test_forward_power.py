from power.lowz_power import compute_power_forward
import astropy.table as tbl
import numpy as np

forward_file ='/home/vikram/output_power/test/danforth_forward_models_min_SN_10_Gamma_0.09000_Nran_010000_seed_42.fits'

Power_table = tbl.Table( compute_power_forward(forward_data=forward_file,  use_metalmasking=True),
                         names=['LS K-modes', 'LS Power(noise and wind corr)', 'N_used', 'C', 'z'])
Power_table.add_column(tbl.Column(np.sqrt(np.diag(Power_table['C'])), name='sigma'))
out = Power_table['z', 'LS K-modes', 'LS Power(noise and wind corr)', 'sigma', 'N_used']

oufile_name = '/home/vikram/output_power/test/power_forward_tng_Gamma_0.09000_Nran_010000_seed_42.txt'
out.write(oufile_name, format="ascii.commented_header", overwrite=True)