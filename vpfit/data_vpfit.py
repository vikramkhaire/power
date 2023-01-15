import matplotlib as mpl
mpl.use ('Agg')

import numpy as np
import astropy.table as tab
from pypowerspec.data_readin import data_readin

from astropy.io import ascii
import os


def just_write_spectra(data, res_array, outpath, outfilename, select=None):

    # setting names for output dir
    outpathname = outpath + '/gamma_' + outfilename+ '_res_cos'

    print('stroring in dir: ', outpathname)
    if not os.path.isdir(outpathname):
        os.mkdir(outpathname)
        print('creating dir: ', outpathname)

    # continuum and mask
    number_of_pixels = len(data['Wave'][0])
    continuum = np.ones(number_of_pixels)
    mask = np.zeros(number_of_pixels).astype(int) # actual masks are used to remove contaminations

    data_length = len(data)
    for i in range(data_length):
        ndata = tab.Table([data['Wave'][i], data['Flux'][i], data['Noise'][i], continuum, mask])
        # removing Nan valuse
        ndata = ndata[~np.isnan(ndata['col1'])]

        ndata.meta['comments'] = ['#']  # just this
        ascii.write(ndata, outpathname + '/spectrum_{}_COS{}.dat'.format(data['Obj'][i], res_array[i].lower()),
            format='no_header', overwrite=True,
            comment='')
        print('writing', i, data['Obj'][i])
        # Hector's code need lowercase lp
        if not i % 10:
            print('writing {} th spectrum in {}'.format(i, outpathname))

    return


def write_each_bin(datadir, outpath, min_z= 0.005, max_z= 0.06, minsn = -0.1, cutrange=[1050, 1180], minrange = 1e-10):

    data, resolution_array, z_qso_array = data_readin(use_metalmasking=True, min_z=min_z, max_z=max_z,
        dataset='COS_data', cutrange = cutrange, minrange = minrange,
        minsn=minsn, path=datadir, fill_with_noise=True, use_emissionmasking=True)

    # store data (which has original mask column)
    file= 'data_read_z_{:0.3f}_{:0.3f}.fits'.format(min_z, max_z)
    data.write(outpath+'/'+file, overwrite = True)

    # use_metalmasking and use_emissionmasking needs to be True

    outfilename = 'z_{:0.3f}_{:0.3f}'.format(min_z, max_z)

    print(len(data), ' .... data that is read')

    just_write_spectra(data=data, res_array=resolution_array, outpath=outpath, outfilename= outfilename)


    return


def write_hsla_spectra_in_correct_format(hsla_path, outpath):

    master = tab.Table.read(hsla_path + '/' + 'masterfile.fits')

    for i in range(len(master)):
         qname = master['obj'][i]
         spectrum_filename =  hsla_path + '/' + qname + '_coadd_G130M_final_lpALL_continuum.fits'
         data = tab.Table.read(spectrum_filename)

         new_filename = outpath + '/' + qname + '.dat'
         ndata = tab.Table([data['WAVE'], data['FLUX'], data['ERROR'], data['Conti_spline'], data['mask']])
         # removing Nan valuse
         ndata = ndata[~np.isnan(ndata['FLUX'])]

         ndata.meta['comments'] = ['resvel 15.0']  # just this
         ascii.write(ndata, new_filename, format='no_header', overwrite=True, comment='')
         print('writing', i, 'spectrum', qname)

    # write masterfile
    master_file_new = outpath + '/'  + 'masterfile.txt'
    ascii.write(master, master_file_new, format='no_header', overwrite=True, comment='')

    return

"""
hsla_path = '/mnt/quasar/vikram/hsla/hsla'
outpath = '/mnt/quasar/vikram/hsla/mw_format'

write_hsla_spectra_in_correct_format(hsla_path =  hsla_path, outpath = outpath)

datadir = outpath
vpfit_outpath = '/mnt/quasar/vikram/hsla/vpfit_all'
write_each_bin(datadir=datadir, outpath=vpfit_outpath, min_z= 0.005, max_z= 0.16)
"""

hsla_path = '/home/vikram/output_power/data'
outpath = '/home/vikram/output_power/mw_format'

write_hsla_spectra_in_correct_format(hsla_path =  hsla_path, outpath = outpath)

datadir = outpath +'/'
vpfit_outpath = '/home/vikram/output_power/vpfit_all'
write_each_bin(datadir=datadir, outpath=vpfit_outpath, min_z= 0, max_z= 0.6)


