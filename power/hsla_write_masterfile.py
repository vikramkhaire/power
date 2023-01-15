import astropy.table as tab
import numpy as np
import glob
import math


data_path = '/home/vikram/output_power/data'

d  =  tab.Table.read(data_path + '/' + 'quasar_probing_low_z_lyaf_Jan14.txt', format ='ascii')
qname = d['qname_lyaf']

files =  glob.glob(data_path + '/*.fits')

qname_stored = []
for i in files:
    name  = i.split('/')[-1].split('_coadd')[0]
    qname_stored.append(name)


#--- things to write in the
obj = []
redshift = []
sn = []  # did not store; cut was applied to be > 5
flag = [] # everthing is 1
res = []  # fix 15 km/s
lpin = [] # fixing it to lifetime position 1

quality_array = []
qname =  list(qname)
for j in qname_stored:
    try:
        ind = qname.index(j)


        # calculating SN
        fits_filename = data_path +'/' + j + '_coadd_G130M_final_lpALL_continuum.fits'
        data = tab.Table.read(fits_filename)
        data = data[data['mask']==0]

        count_short = 0
        if len(data) >10 :
            SigN = np.median(data['FLUX'] / data['ERROR']) * np.sqrt(3)
            SigN = round(SigN, 1)

            print(SigN, j)
            sn.append(SigN)

            obj.append(j)
            redshift.append(d['zem_q_lyaf'][ind])

            """
            # --- keeping the const numbers
            sn.append(10)
            # sor flags in hsla data by Sapna
            # flag 1 for SNR < 3 and flag 2 for 3<SNR <5
            data = tab.Table.read(data_path + '/{}_coadd_G130M_final_lpALL_continuum.fits'.format(j))
            quality = (len(data)- len(data[data['GAP_FLAGS']==6.0])) / (len(data[data['GAP_FLAGS']==0.0]))
            #if quality > 3.1: # decided depending upon the histogram of quality to include many qsos for check (original was just 2)
            if quality > 2:  # fiducial
                flag.append(0)
                #print(quality, '->0')
            else:
                flag.append(1)
                #print(quality, '->1')
            """
            res.append(15)
            lpin.append('LP1')
            quality_array.append(1)
            flag.append(1)
        else:
            count_short +=1
            print('bad spctra #', count_short, j)


    except:
        print(j, 'not present in the file quasar_probing_low_z_lyaf.txt')

# store masterfile
master_table = tab.Table([obj, redshift, sn, flag, res, lpin], names = ('obj', 'z', 'sn', 'flag', 'res', 'lp'))
#print(master_table)

masterfile_name = data_path + '/masterfile.fits'
master_table.write(masterfile_name, overwrite = True)
