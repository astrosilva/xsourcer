# initialize ciao env
import os
from math import trunc
import numpy as np
import matplotlib.pylab as plt
from astropy.convolution import Box1DKernel, convolve
from astropy.io import fits

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def met_jd(met):
    met_mjd = (met / 86400) + 50814.0  # Terrestrial time (ref 1998.0)
    mjd_jd = met_mjd + 2400000.5       # ref January 1, 4713 BCE
    jd_utc = mjd_jd - (32.184 + 32) / 86400  # TT-UTC + TAI-UTC

    return jd_utc

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#target = os.getcwd().rsplit('/')[-1]
#obsid = '3392_comp'
target = os.getcwd().rsplit('/')[-2]

os.chdir('multi_obs/')

bands = [('total', '0.1-10'),('soft', '0.1-2'), ('hard', '2-10')]
energies = [b + ' keV' for (i,b) in bands]

fig, ax = plt.subplots(3, 1, figsize=(12,5), sharex=True)
plt.subplots_adjust(hspace=0.0)
for idx, (band, energy) in enumerate(zip(bands, energies)):
    hdul = fits.open(f'lc_{band[0]}_{band[1]}.fits')[1]
    counts = hdul.data['net_rate']
    error = hdul.data['err_rate']
    
    tstart = met_jd(hdul.header['TSTART'])
    tstop = met_jd(hdul.header['TSTOP'])
    dt = np.linspace(tstart, tstop, len(counts))

    ax[idx].errorbar(dt, counts, yerr=error, color='black', marker='None', label=energy, lw=0.8)
    ax[idx].legend(loc='upper right')
    ax[idx].set_xlim(tstart, tstop)

    if idx == 1:
        ax[idx].set_ylabel(r'COUNT_RATE   [s$^{-1}$]', labelpad=10)
    if idx == 2:
        ax[idx].set_xlabel('TIME [BJD]')
        ax[idx].locator_params(axis='x', nbins=5)

        
def phase(tstart, t0, period):
    return (tstart - t0) / period

#ref_phase = trunc(phase(tstart, T0, P))
#ref_phase = np.floor(phase(ts, T0, P)).astype(int)
#phase_ticks = np.array([0.8, 1.0, 1.2, 1.4, 1.6])
#phase_ticks = phase(time_ticks, T0, P) - ref_phase
#time_ticks = (phase_ticks * P) + T0 + (ref_phase * P)

#ax_top = ax[0].twiny()
#ax_top.set_xlim(ax[0].get_xlim())
#ax_top.set_xticks(time_ticks)
#ax_top.set_xticklabels([f'{p:.2f}' for p in phase_ticks])
#ax_top.set_xlabel('ORBITAL PHASE', labelpad=10)

# Eclipse phase range
#eclipse_start_phase = 0.95
#eclipse_end_phase = 0.05  # wraps around

# Eclipse in time — generate all eclipse windows covering the time range
#n_cycles = int(np.ceil((TSTOP - TSTART) / P)) + 1
#for i in range(n_cycles):
#    cycle_start = T0 + (ref_phase + i) * P
#    eclipse_start_time = cycle_start + eclipse_start_phase * P
#    eclipse_end_time = cycle_start + (1 + eclipse_end_phase) * P  # add 1 to wrap

    # Only plot if eclipse overlaps with the time window
#    if eclipse_end_time >= TSTART and eclipse_start_time <= TSTOP:
#        for a in ax:
#            a.axvspan(eclipse_start_time, eclipse_end_time, color='red', alpha=0.2)

#plt.savefig('lightcurve.pdf')
plt.show()
exit()


# Define spin parameters
SPIN_PERIOD = 0.014517757
T0_SPIN = 2456977.1057

one = read_file('FO_Aqr_0.3-2_lc.fits')
two = read_file('FO_Aqr_2-4_lc.fits')
three = read_file('FO_Aqr_4-6_lc.fits')
four = read_file('FO_Aqr_6-7_lc.fits')
bands = [one, two, three, four]
ranges = ['0.3-2 keV','2-4 keV', '4-6 keV', '6-7 keV']
fig, ax = plt.subplots(4, 1, figsize=(10, 6), sharex=True)
plt.subplots_adjust(hspace=0)

for idx, (band, range) in enumerate(zip(bands, ranges)):

    rates = band.get_column('net_rate').values
    times = np.linspace(TSTART, TSTOP, len(rates))
    
    spin_phase = ((times - T0_SPIN) / SPIN_PERIOD) % 1
    # Duplicate for two cycles
    full_phase = np.concatenate([spin_phase, spin_phase + 1])
    full_rate = np.concatenate([rates, rates])
    
    sorted_idx = np.argsort(full_phase)
    full_phase = full_phase[sorted_idx]
    full_rate = full_rate[sorted_idx]
    box_kernel = Box1DKernel(15)
    full_rate = convolve(full_rate, box_kernel)

    ax[idx].plot(full_phase, full_rate, lw=0.8, color='black', label=range)
    ax[idx].legend(loc='upper right')
    ax[0].set_xlim(0, 2)
    ax[0].set_ylim(0, 2)
    if idx > 0:
        ax[idx].set_ylim(0, 1.99)

fig.text(0.055, 0.5, r'COUNT_RATE   [s$^{-1}$]', va='center', rotation='vertical', fontsize=11)
ax[-1].set_xlabel('SPIN PHASE', labelpad=10, fontsize=11)
plt.savefig('spin_phase.pdf')
plt.show()
