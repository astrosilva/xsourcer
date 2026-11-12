import os
import numpy as np
from astropy.io import fits
from scipy.signal import lombscargle
import matplotlib.pyplot as plt
from pycrates import read_file
from ciao_contrib.runtool import *
from subprocess import Popen

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def lomb_scargle(times, rates):
    rates -= np.mean(rates)    
    frequencies = np.arange(fmin, fmax, 1e-6)
    angular = 2 * np.pi * frequencies

    power = np.zeros_like(angular)
    for i, w in enumerate(angular):
        tau = (1 / (2 * w)) * np.arctan(np.sum(np.sin(2 * w * times)) / np.sum(np.cos(2 * w * times)))

        num_t1 = np.sum(rates * np.cos(w * (times - tau)))**2
        den_t1 = np.sum(np.cos(w * (ttime - tau))**2)
        num_t2 = np.sum(rates * np.sin(w * (times - tau)))**2
        den_t2 = np.sum(np.sin(w * (ttime - tau))**2)
        power[i] = ((num_t1 / den_t1) + (num_t2 / den_t2)) / 2
    
    return frequencies, power


def rayleigh(times):
    N = len(times)
    frequencies = np.arange(fmin, fmax, 1e-6)
    
    power = np.zeros_like(frequencies)
    for i, f in enumerate(frequencies):
        phi = 2 * np.pi * f * times
        power[i] = (2 / N) * (np.sum(np.cos(phi))**2 + np.sum(np.sin(phi))**2)

    return frequencies, power

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

target = 'CXOGC_J174543.4-285841'
#obsid = '3392_comp'
#os.chdir(f'{obsid}/reduction/')

os.chdir(f'{target}/multi_obs/')

#total = read_file('lc_total_0.1-10.fits')

#trate = total.get_column('net_rate').values
#ttime = total.get_column('dt').values

fmin, fmax = 3e-5, 1e-2
#lfreq, lpow = lomb_scargle(ttime, trate)

ztime = read_file('merged_evt.fits[sky=region(source.reg)]').get_column('TIME').values
zfreq, zpow = rayleigh(ztime)

#lpeak = lfreq[np.where(lpow == max(lpow))][0]
zpeak =	zfreq[np.where(zpow == max(zpow))][0]
#llpeak = 1 / lpeak
zzpeak = 1 / zpeak

# Plot
fig, ax = plt.subplots(2, 1, figsize=(12, 4), sharex=True)
fig.subplots_adjust(hspace=0)

#ax[0].plot(lfreq, lpow, lw=0.9, color='black', label='Lomb-Scargle')
#ax[0].axvline(lpeak, ls='dashed', lw=0.9, color='red', label=f'{llpeak:.2f} s')
ax[0].set_xlim(fmin, fmax)
ax[0].set_ylim(bottom=0)
ax[0].set_xscale('log')
ax[0].set_ylabel('Power', fontsize=11, labelpad=10)
ax[0].set_title(f'Lomb-Scargle Periodogram : {target}', fontsize=13)
#ax[0].legend()

ax[1].plot(zfreq, zpow, lw=0.9, color='black', label=r'$Z^2_1$ (Rayleigh)')
ax[1].axvline(zpeak, ls='dashed', lw=0.9, color='red', label=f'{zzpeak:.2f} s')
ax[1].set_ylim(bottom=0)
ax[1].set_xlabel('Frequency [Hz]', fontsize=11, labelpad=5)
ax[1].set_ylabel('Power', fontsize=11, labelpad=10)
ax[1].legend()

#plt.savefig('periodogram.pdf')
plt.show()











# PHASE FOLDED LIGHTCURVE
#evt2 = str([i for i in os.listdir() if i.endswith('bary_evt2.fits')][0])
#start = fits.open(evt2)[1].header['TSTART']

#folder = ['pfold',
#          f'infile={evt2}',
#          'outfile=flc_total_0.1-10.fits',
#          'periodgrid=5200:5400:0.01',
#          f'tzero={start}',
#          'clobber=yes']
#verbose(folder)
#exit()

#folded = read_file('flc_total_0.1-10.fits')
#period = folded.get_column('period').values
#mean_rate = folded.get_column('mean_rate').values
#sigma_rate = folded.get_column('sigma_rate').values

#fig, ax = plt.subplots(2, 1, figsize=(12, 4), sharex=True)
#plt.subplots_adjust(hspace=0.0)
#ax[0].plot(period,mean_rate, marker="None")
#ax[0].set_ylabel("Mean Count Rate (cts/s)")

#ax[1].plot(period,sigma_rate, marker="None")
#ax[1].set_ylabel("Stddev. Count Rate (cts/s)")
#ax[1].set_xlabel("Period (s)")

#print(period[np.where(mean_rate == max(mean_rate))], period[np.where(sigma_rate == max(sigma_rate))])
#plt.show()
#period = 5328.66#.288#input('Period: ')

#folded = read_file(f'flc_total_0.1-10.fits[period={period}]')
#x=folded.get_column("phase_mid").values.flatten()
#y=folded.get_column("rate").values.flatten()

#plt.plot(x,y,marker="None", drawstyle="steps-mid")
#plt.xlabel("Phase")
#plt.ylabel("Count Rate (cts/s)")
#plt.show()


