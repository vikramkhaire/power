import numpy as np
from matplotlib import pyplot
from astropy.stats import LombScargle
from scipy import interpolate
import astropy.table as tbl
import astropy.units as u
import linetools.spectra.lsf as lt_lsf
import astropy.constants as const
from numpy.fft import rfft,rfftfreq


# setting up the random seeds
np.random.seed(0)

#Useful Constants
Ly_alpha = 1215.67  # Angstroms
c = const.c.to('km/s').value


def gaussian_lsf(k,R):
    return (np.exp(-((k*R)**2)/2))

def stis_lsf(k,wave):
    inst = 'STIS'
    slit = '0.2x0.2'
    grating = 'E230M'
    cen = 2269
    l = lt_lsf.LSF({'name': inst, 'grating': grating, 'cen_wave': str(cen), 'slit': slit})
    wave_even = np.linspace(np.min(wave), np.max(wave), 2 * len(wave)) * u.angstrom
    kern = l.get_lsf(wave_even, 'cubic')
    Wlsf = np.abs(rfft(kern))
    klsf = rfftfreq(len(wave_even), 299792.458 * np.diff(wave_even)[0] / (np.mean(wave_even))) * np.pi * 2
    Wint = interpolate.interp1d(klsf, Wlsf)
    return (Wint(k))

  
def cos_lsf(k,wave,grating,LP):
    inst = 'COS'
    if grating=='G130M':
        cen = 1300
    elif grating=='G160M':
        cen = 1600
    l = lt_lsf.LSF({'name': inst, 'grating': grating, 'cen_wave': str(cen), 'life_position': LP})
    wave_even = np.linspace(np.min(wave), np.max(wave), 2 * len(wave)) * u.angstrom  # why?
    kern = l.get_lsf(wave_even, 'cubic')
    Wlsf = np.abs(rfft(kern))
    klsf = rfftfreq(len(wave_even), 3e5 * np.diff(wave_even)[0] / (np.mean(wave_even))) * np.pi * 2  # why?
    Wint = interpolate.interp1d(klsf, Wlsf)
    return (Wint(k))
  

def window(k_in,delv,lsf=gaussian_lsf,**kwd):
    k = np.asarray(k_in)
    exp_term = lsf(k,**kwd)
    sin_term = (np.sin(k*delv/2))/(k*delv/2)
    sin_term[k==0]=1
    return sin_term*exp_term
  
def average(bins,power,k):
    k_out = []
    power_out = []
    N_use_out = []
    for minbin, maxbin in zip(bins[:-1],bins[1:]):
        filt = (k <= maxbin) & (k> minbin)
        N_use=np.sum(filt)
        N_use_out.append(N_use)
        k_use = np.mean(k[filt])
        power_use = np.mean(power[filt])
        k_out.append(k_use)
        power_out.append(power_use)
    return np.array(k_out), np.array(power_out), np.array(N_use_out)
    
###########################

def compute_power(use_metalmasking=False,wavelim=[1050,1180],zbin=[0.005,0.171875],no_lsf_correction=False,
                  fill_with_noise=False, only_130M=False, only_160M=False, use_milkyway_metals=False, use_igm_metals=False):
    
    # for central_z in redshift_bins:
    Lomb_k = []
    Lomb_power = []
    z_abs_all = []
    path = '/home/vikram/Dropbox/power_spectrum_idl_cleaned_up/'
    # Looping through QSO spectra #
    for obj, zqso, res, lp in zip(obj_name, redshift, Resolution, LP):
        data = np.loadtxt(path + 'danforth_data_MW_format_new/modified_MW_format/{}.dat'.format(obj.astype(str)),
                          comments='resvel')
        # print (obj, zqso)


        ## Reading Data ##
        wavelength = data[:, 0]
        flux = data[:, 1]
        sigma_F = data[:, 2]
        continuum = data[:, 3]
        mask = data[:, 4]

        # cutting data into rest frame range 1050 Ang to 1180Ang #
        cut2 = np.logical_and(wavelength > wavelim[0] * (1 + zqso), wavelength <= wavelim[1] * (1 + zqso))
        # cutting data in redshift range central_z -0.1 and central_z +0.1 #
        cut3 = np.logical_and(wavelength > (1 + zbin[0]) * Ly_alpha, wavelength <= (1 + zbin[1]) * Ly_alpha)
        # combining cuts with masking. here we include for mask==5 data#
        cut = cut2 & cut3 & (((mask==0) | (mask==5) | (mask==7)) if use_metalmasking==False else mask==0)
        if use_milkyway_metals==True:
            cut = cut2 & cut3 & ((mask==0) | (mask==7))
        if use_igm_metals==True:
            cut = cut2 & cut3 & ((mask==0) | (mask==5))

        ## This only allows mask==0 or 5 into the LS data set ##
        if fill_with_noise:
            cutall0 = (mask==0) | (mask==5) | (mask==6) | (mask==7)
            cutall = cutall0 & cut2 & cut3
            wavelength = wavelength[cutall]
            flux = flux[cutall]
            sigma_F = sigma_F[cutall]
            continuum = continuum[cutall]
            mask = mask[cutall]
            remove = (mask==5) | (mask==6)
            flux[remove] = continuum[remove] + np.random.randn(np.sum(remove)) * sigma_F[remove]
        else:
            wavelength = wavelength[cut]
            flux = flux[cut]
            sigma_F = sigma_F[cut]
            continuum = continuum[cut]
            mask = mask[cut]

        # Remove any spectra 'too short' (i.e. less than mentioned wavelength)#
        minrange = (np.log(wavelim[1]) - np.log(1050.)) / 10.
        if (len(wavelength)==0):
            continue
        elif (np.max(np.log(wavelength)) - np.min(np.log(wavelength)) < minrange):
            continue

        ## Gathering/Converting data points into km/s ##
        z_absorber = (wavelength / Ly_alpha) - 1
        z_bar = np.mean(z_absorber)
        F_bar = np.mean(flux / continuum)
        velocity_data = c * np.log(wavelength / Ly_alpha / (1 + zqso))

        tauevo = -np.log(F_bar)  # -np.log(0.963)#tau0 * ((1+z_absorber)/(1+z0))**beta + C  #adjusted optical depth
        dF_list = flux / (np.exp(-tauevo)) / continuum - 1
        normalized_sigma = sigma_F / (np.exp(-tauevo)) / continuum

        print('{} has a redshift of {} and Fbar of {:.4f} with zbar of {:.4f}'.format(obj, zqso, F_bar, z_bar))
        # Length of Spectra in velocity units #
        resfactor = np.median(velocity_data[1:] - velocity_data[:-1])

        #### Unmasked LS ###
        LS = LombScargle(velocity_data, dF_list, fit_mean=False, center_data=True)
        LS_k, LS_raw_power = LS.autopower(method='fast', normalization='psd', samples_per_peak=1, nyquist_factor=1,
                                          minimum_frequency=0)
        LS_k = 2.0 * np.pi * LS_k[1:]
        LS_raw_power = LS_raw_power[1:]
        LS_norm_raw = resfactor * LS_raw_power  # /(len(dF_list))
        Lomb_k.append(LS_k)
        ## Generate Noise from normalized_sigma ##
        LS_noisy_power = []
        for jjj in range(100):
            noise = (normalized_sigma) * np.random.standard_normal(len(dF_list))
            LS_noise = LombScargle(velocity_data, noise, fit_mean=False, center_data=True)
            LS_noise_k, LS_noise_power = LS_noise.autopower(method='fast', normalization='psd', samples_per_peak=1,
                                                            nyquist_factor=1, minimum_frequency=0)
            LS_noise_k = 2.0 * np.pi * LS_noise_k[1:]
            LS_noise_power = LS_noise_power[1:]
            LS_norm_noise = resfactor * LS_noise_power  # /len(noise)
            # print (LS_noise_k - LS_k)
            LS_noisy_power.append(LS_norm_noise)
        LS_noise_mean = np.mean(LS_noisy_power, axis=0)
        Raw_diff_Noise = LS_norm_raw - LS_noise_mean
        if not no_lsf_correction:
            ith_power = Raw_diff_Noise / (window(LS_k, resfactor, lsf=cos_lsf, wave=wavelength, LP=lp,
                                                 grating=('G130M' if np.median(wavelength) < 1460.0 else 'G160M')) ** 2)
            if only_130M:
                ith_power = Raw_diff_Noise / (window(LS_k, resfactor, lsf=cos_lsf, wave=wavelength, LP=lp, grating=(
                    'G130M' if np.median(wavelength) < 1460.0 else 'G160M')) ** 2)
                print('G130M' if np.median(wavelength) < 1500.0 else 'G160M')
            if only_160M:
                ith_power = Raw_diff_Noise / (
                            window(LS_k, resfactor, lsf=cos_lsf, wave=wavelength, LP=lp, grating=('G160M')) ** 2)
        else:
            ith_power = Raw_diff_Noise
        Lomb_power.append(ith_power)
        z_abs_all.append(z_absorber)
    try:
        LS_k = np.concatenate(np.array(Lomb_k))
    except ValueError:
        return [[], [], [], [], []]
    LS_power = np.concatenate(Lomb_power)

    bins = np.logspace(np.log10(0.00012596187096205269 / np.sqrt(10)), np.log10(1.25962), 46)
    # FFT_masked_k_mean, FFT_masked_power_mean = average(bins/(2.0*np.pi), power, k)

    LS_k_mean, LS_power_mean, N_used = average(bins, LS_power, LS_k)
    ### We can randomly sample N spectra from the power spectra stored above ##
    bootstrap_sample = []
    for j in range(1000):
        randind = np.random.randint(0, len(Lomb_power), size=len(Lomb_power))
        randPower = np.concatenate([Lomb_power[r] for r in randind])
        randk = np.concatenate([Lomb_k[r] for r in randind])
        rand_k_mean, rand_power_mean, N_used_boot = average(bins, randPower, randk)
        bootstrap = np.outer(LS_k_mean * (rand_power_mean - LS_power_mean) / np.pi,
                             LS_k_mean * (rand_power_mean - LS_power_mean) / np.pi)
        bootstrap_sample.append(bootstrap)
    Cij = np.nanmean(bootstrap_sample, axis=0)
    return [LS_k_mean, LS_k_mean * LS_power_mean / np.pi, N_used, Cij,
            [np.mean(np.concatenate(z_abs_all))] * len(LS_k_mean)]




path = '/home/vikram/Dropbox/power_spectrum_idl_cleaned_up/'
##### Reading in COS QSO and z ####
obj_name, redshift, sn, flag, res, lpin = zip(
    *np.genfromtxt(path + 'danforth_data_MW_format_new/masterfile.txt', comments='#',
                   dtype=['<U30', np.float, np.float, np.int, np.float, 'U3']))

obj_name, redshift, sn, flag, res, lpin = np.array(obj_name), np.array(redshift), np.array(sn), np.array(
    flag), np.array(res), np.array(lpin)

###### Pick Out COS Sample (i.e. flags) and Put in RES #########

obj_name = obj_name
redshift = redshift
Resolution = res / 2.35  # Why the division?

# obj_name = data2[:,0]
filtered = (sn >= 10) & (flag==1)
redshift = redshift[filtered]
obj_name = obj_name[filtered]
sn = sn[filtered]
flag = flag[filtered]
Resolution = Resolution[filtered]
LP = [l[-1] for l in lpin[filtered]]
# LP=[l.decode()[-1] for l in lpin]


zbinarr=[[0.005, 0.06], [0.06, 0.16], [0.16, 0.26], [0.26, 0.36], [0.36, 0.48]]
#zbinarr=[[0.16, 0.26]]

#zbinarr=[[0.36, 0.48]]
#zbinarr=[[0.005, 0.26]]
#zbinarr=[[0.001, 0.06], [0.06, 0.16], [0.16, 0.26], [0.26, 0.36],[0.36, 0.48]]
#zbinarr=[[0.005, 0.04], [0.005, 0.05], [0.005, 0.06], [0.005, 0.07], [0.005, 0.08], [0.005, 0.09], [0.005, 0.1], [0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.01, 0.11], [0.11, 0.21], [0.21, 0.31], [0.31, 0.41], [0.02, 0.12], [0.12, 0.22], [0.22, 0.32], [0.32, 0.42], [0.03, 0.13], [0.13, 0.23], [0.23, 0.33], [0.33, 0.43], [0.04, 0.14], [0.14, 0.24], [0.24, 0.34], [0.34, 0.44], [0.05, 0.15], [0.15, 0.25], [0.25, 0.35], [0.35, 0.45], [0.06, 0.16], [0.16, 0.26], [0.26, 0.36], [0.36, 0.46], [0.07, 0.17], [0.17, 0.27], [0.27, 0.37], [0.37, 0.47], [0.08, 0.18], [0.18, 0.28], [0.28, 0.38], [0.38, 0.48], [0.09, 0.19], [0.19, 0.29], [0.29, 0.39], [0.39, 0.49]] 
waverangearr=[[1050, 1180]]
#[[930,1195]]
#,[950,1180],[0,1180]]
for zbin in zbinarr:
  for waverange in waverangearr:
      outpath = '/home/vikram/output_power/'
      Power_table = tbl.Table(compute_power(use_metalmasking=True,zbin=zbin,wavelim=waverange),names=['LS K-modes','LS Power(noise and wind corr)','N_used', 'C','z'])
      Power_table.add_column(tbl.Column(np.sqrt(np.diag(Power_table['C'])),name='sigma'))
      out=Power_table['z','LS K-modes','LS Power(noise and wind corr)','sigma','N_used']
      out.write(outpath+'z0Power_metals_masked_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)
    

  #  Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,only_160M=True, zbin=zbin,wavelim=waverange),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
  #  Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
  #  out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
  #  out_masked.write('Joint_Output/test_grating/160Mz0Power_metals_masked_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)
    

    #Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,zbin=zbin,wavelim=waverange,no_lsf_correction=True),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
    #Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
    #out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
    #out_masked.write('Joint_Output/z0Power_metals_masked_no_lsf_correct_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)
  
    #Power_table_masked = tbl.Table(compute_power(use_metalmasking=True,zbin=zbin,wavelim=waverange,fill_with_noise=True),names=['LS K-modes','masked LS Power (noise and wind corr)','N_used','C','z'])
    #Power_table_masked.add_column(tbl.Column(np.sqrt(np.diag(Power_table_masked['C'])),name='sigma'))
    #out_masked=Power_table_masked['z','LS K-modes','masked LS Power (noise and wind corr)','sigma','N_used']
    #out_masked.write('test_filled_masks/z0Power_filled_masks_z{:.2f}-{:.2f}_wave{}-{}.txt'.format(zbin[0],zbin[1],waverange[0],waverange[1]),format="ascii.commented_header", overwrite=True)

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


#pyplot.figure()
#pyplot.loglog(2.0*np.pi*LS_masked_k_mean, 2.0*np.pi*LS_masked_k_mean*LS_masked_power_mean/np.pi,'r-', label='LS (masking 5, noise and wind corr)')
#pyplot.xlabel('k s/km')
#pyplot.ylabel('kP(k)/pi')


#pyplot.legend()
#pyplot.show()

