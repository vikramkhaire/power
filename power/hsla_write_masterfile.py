import astropy.table as tab
import numpy as np
import glob


data_path = '/home/vikram/output_power/hsla'

d  =  tab.Table.read(data_path + '/' + 'quasar_probing_low_z_lyaf.txt', format ='ascii')
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

qname =  list(qname)
for j in qname_stored:
    try:
        ind = qname.index(j)
        obj.append(j)
        redshift.append(d['zem_q_lyaf'][ind])
        # --- keeping the const numbers
        sn.append(10)
        # sor flags in hsla data by Sapna
        # flag 1 for SNR < 3 and flag 2 for 3<SNR <5
        data = tab.Table.read(data_path + '/{}_coadd_G130M_final_lpALL_continuum.fits'.format(j))
        if len(data)/len(data[data['GAP_FLAGS']==0.0]) > 2:
            flag.append(0)
        else:
            flag.append(1)
        res.append(15)
        lpin.append('LP1')
    except:
        print(j, 'not present in the file quasar_probing_low_z_lyaf.txt')

# store masterfile
master_table = tab.Table([obj, redshift, sn, flag, res, lpin], names = ('obj', 'z', 'sn', 'flag', 'res', 'lp'))
print(master_table)

masterfile_name = data_path + '/masterfile.fits'
master_table.write(masterfile_name, overwrite = True)
