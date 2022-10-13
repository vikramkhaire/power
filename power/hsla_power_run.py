from power.lowz_power import compute_power_hsla
import astropy.table as tbl
import numpy as np

data_path = '/home/vikram/output_power/hsla'

zbinarr = [[0.005, 0.06], [0.06, 0.16], [0.16, 0.26], [0.26, 0.36]]
#zbinarr = [ [0.06, 0.16]]

waverangearr = [[1050, 1180]]

for zbin in zbinarr:
    print('working in z=', zbin)
    for waverange in waverangearr:
        outpath = '/home/vikram/output_power/hsla_power/'
        Power_table = tbl.Table(
            compute_power_hsla(data_path=data_path, use_metalmasking=True, zbin=zbin, wavelim=waverange, fill_with_noise = False),
            names=['LS K-modes', 'LS Power(noise and wind corr)', 'N_used', 'C', 'z'])
        Power_table.add_column(tbl.Column(np.sqrt(np.diag(Power_table['C'])), name='sigma'))
        out = Power_table['z', 'LS K-modes', 'LS Power(noise and wind corr)', 'sigma', 'N_used']
        out.write(outpath + 'Power_metals_masked_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0], zbin[1], waverange[0],
                                                                                        waverange[1]),
                  format="ascii.commented_header", overwrite=True)

