import sys
import matplotlib as mpl
mpl.use('Agg')
import os
import vpfit_controller as vpc

path = '/mnt/quasar/vikram/hsla'
os.chdir(path)
vpc.vpfit_godzilla.vpfit_main(0, -1, 1, 0, -1, data = False, custom_lsf = True)