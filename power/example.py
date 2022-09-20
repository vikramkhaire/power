from power.lowz_power import compute_power
import astropy.table as tbl

masterfile = '/home/vikram/Dropbox/power_spectrum_idl_cleaned_up/danforth_data_MW_format_new/masterfile.txt'
data_path = '/home/vikram/Dropbox/power_spectrum_idl_cleaned_up/danforth_data_MW_format_new/modified_MW_format'

zbinarr = [[0.005, 0.06], [0.06, 0.16], [0.16, 0.26], [0.26, 0.36], [0.36, 0.48]]
# zbinarr=[[0.16, 0.26]]

# zbinarr=[[0.36, 0.48]]
# zbinarr=[[0.005, 0.26]]
# zbinarr=[[0.001, 0.06], [0.06, 0.16], [0.16, 0.26], [0.26, 0.36],[0.36, 0.48]]
# zbinarr=[[0.005, 0.04], [0.005, 0.05], [0.005, 0.06], [0.005, 0.07], [0.005, 0.08], [0.005, 0.09], [0.005, 0.1], [0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.01, 0.11], [0.11, 0.21], [0.21, 0.31], [0.31, 0.41], [0.02, 0.12], [0.12, 0.22], [0.22, 0.32], [0.32, 0.42], [0.03, 0.13], [0.13, 0.23], [0.23, 0.33], [0.33, 0.43], [0.04, 0.14], [0.14, 0.24], [0.24, 0.34], [0.34, 0.44], [0.05, 0.15], [0.15, 0.25], [0.25, 0.35], [0.35, 0.45], [0.06, 0.16], [0.16, 0.26], [0.26, 0.36], [0.36, 0.46], [0.07, 0.17], [0.17, 0.27], [0.27, 0.37], [0.37, 0.47], [0.08, 0.18], [0.18, 0.28], [0.28, 0.38], [0.38, 0.48], [0.09, 0.19], [0.19, 0.29], [0.29, 0.39], [0.39, 0.49]]
waverangearr = [[1050, 1180]]
# [[930,1195]]
# ,[950,1180],[0,1180]]
for zbin in zbinarr:
    for waverange in waverangearr:
        outpath = '/home/vikram/output_power/'
        Power_table = tbl.Table(
            compute_power(masterfile=masterfile, data_path=data_path, use_metalmasking=True, zbin=zbin,
                          wavelim=waverange),
            names=['LS K-modes', 'LS Power(noise and wind corr)', 'N_used', 'C', 'z'])
        Power_table.add_column(tbl.Column(np.sqrt(np.diag(Power_table['C'])), name='sigma'))
        out = Power_table['z', 'LS K-modes', 'LS Power(noise and wind corr)', 'sigma', 'N_used']
        out.write(outpath + 'z0Power_metals_masked_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0], zbin[1], waverange[0],
                                                                                        waverange[1]),
                  format="ascii.commented_header", overwrite=True)


    #  Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,only_160M=True, zbin=zbin,wavelim=waverange),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
    #  Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
    #  out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
    #  out_masked.write('Joint_Output/test_grating/160Mz0Power_metals_masked_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)


    # Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,zbin=zbin,wavelim=waverange,no_lsf_correction=True),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
    # Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
    # out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
    # out_masked.write('Joint_Output/z0Power_metals_masked_no_lsf_correct_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)

    # Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,zbin=zbin,wavelim=waverange,fill_with_noise=True),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
    # Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
    # out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
    # out_masked.write('test_filled_masks/z0Power_filled_masks_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)

"""
    Power_table_masked = tbl.Table(compute_power(zbin=zbin,wavelim=waverange, use_milkyway_metals=True),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
    Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
    out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
    out_masked.write(outpath+ 'selective_metals/z0Power_with_milkyway_metals_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)

    Power_table_masked = tbl.Table(compute_power(zbin=zbin,wavelim=waverange, use_igm_metals=True ),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
    Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
    out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
    out_masked.write(outpath+'selective_metals/z0Power_with_igm_metals_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)

    Power_table_masked = tbl.Table(compute_power(use_metalmasking=False, zbin=zbin,wavelim=waverange ),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
    Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
    out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
    out_masked.write(outpath+'selective_metals/z0Power_metals_not_masked_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)

"""

#    Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,zbin=zbin,wavelim=waverange, only_130M=True),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
#    Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
#    out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
#    out_masked.write('check_overlap/z0Power_metals_masked_only130M_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)

#    Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,zbin=zbin,wavelim=waverange, only_160M=True),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
#    Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
#    out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
#    out_masked.write('check_overlap/z0Power_metals_masked_only160M_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)

#    Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,zbin=zbin,wavelim=waverange),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
#    Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
#    out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
#    out_masked.write('check_overlap/z0Power_metals_masked_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)


# pyplot.figure()
# pyplot.loglog(2.0*np.pi*LS_masked_k_mean, 2.0*np.pi*LS_masked_k_mean*LS_masked_power_mean/np.pi,'r-', label='LS (masking 5, noise and wind corr)')
# pyplot.xlabel('k s/km')
# pyplot.ylabel('kP(k)/pi')


# pyplot.legend()
# pyplot.show()

