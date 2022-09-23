import astropy.table as tab
import numpy as np
import glob


data_path = ''

d  =  tab.Table.read(data_path + 'quasar_probing_low_z_lyaf.txt', format ='ascii')
qname = d['qname_lyaf']

files =  glob.glob(data_path + '/*.fits')

qname_stored = []
for i in files:
    name  = i.split('_coadd')[0]
    qname_stored.append(name)


#--- things to write in the
obj = []
redshift = []
sn = []  # did not store; cut was applied to be > 5
flag = [] # everthing is zero
res = []  # fix 15 km/s
lpin = [] # fixing it to lifetime position 1

qname =  list(qname)
for j in qname_stored:
    try:
        ind = qname.index(j)
        obj.append(j)
        redshift.append(d['zem_q_lyaf'])
        # --- keeping the const numbers
        sn.append(10)
        flag.append(0)
        res.append(15)
        lpin.append('LP1')
    except:
        print(j, 'not present in the file quasar_probing_low_z_lyaf.txt')

# store masterfile
