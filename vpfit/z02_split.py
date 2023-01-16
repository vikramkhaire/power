
import os
import matplotlib as mpl
mpl.use("Agg")
from enigma.whim.z02_codes.z02_vpfit_forward_models import split
import numpy as np
import glob
os.nice(10)

def run_code(path):
    split(path)
    return

path = '/mnt/quasar/vikram/hsla'
folder = glob.glob(path+ '/gamma*')[0]

split(folder)