# FUNCTIONALS
import logging
from subprocess import run
import numpy as np
import gzip, os, shutil
import matplotlib.pyplot as plt
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
# ASTROPACKS
from astropy.io import fits
from astropy.table import unique, Table
from scipy.signal import find_peaks
from pycrates import read_file
from ciao_contrib.runtool import *
from PyAstronomy.pyasl import coordsDegToSexa as DtS
from source import Source

"""
Process and stack a set of Chandra ACIS observations to produce
final data products using the CIAO software package and procedures
for the purpose of identifying unknown X-ray sources
"""
class Target(Source):
    def query(self, mode=None):
        """
        Creates target directory, queries the Chandra Source Catalog, returns a 
        table of obsids, and reads target RA and Dec
        ----------
        mode : None, 'q'
            Switch for running CSC query. If 'q', only
            target RA and Dec are extracted
        """
        if os.path.exists(f'{self.target_path}/') == True:
            os.chdir(f'{self.target_path}/')

        elif os.path.exists(f'{self.target_path}/') == False:
            os.mkdir(f'{self.target_path}/')
            os.chdir(f'{self.target_path}/')
            self.root = os.getcwd()
            logging.info(f'Created target directory:\n{self.root}\n')

        if mode == None:
            logging.info('Executing query...\n')
            search = ['search_csc',
                      f'pos={self.target}',
                      'radius=0.5',
                      'outfile=query.dat',
                      'bands=broad',
                      'filetype=evt',
                      'catalog=latest',
                      'clobber=yes',
                      'verbose=0']
            self.verbose(search)

            ftable = Table.read('query.dat', format='ascii')
            ftable = ftable[ftable['name'] == f' 2CXO {self.name}']
            ftable = ftable[ftable['instrument'] != ' HRC']
            ftable = ftable['name', 'obsid', 'ra', 'dec']

            ftable.sort('obsid')
            ftable = unique(ftable, keys='obsid', keep='none')
            ftable.write('query.dat', format='ascii', overwrite=True)
            ftable.pprint_all()

        elif mode == 'q':
            ftable = Table.read('query.dat', format='ascii')
        
        ra_deg = ftable['ra'][0]
        dec_deg = ftable['dec'][0]

        sexa = DtS(ra_deg, dec_deg, fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f'))
        self.ra, self.dec = sexa.split(' ')


    def download(self, *obsid_list, mode=None):
        """
        Downloads and gunzips data products for given obsids and
        creates the reduction and repro directories
        ----------
        obsid_list : int, tuple
        mode : None, 'q'
            Switch for downloading data products
        """
        self.obsid_list = obsid_list

        if mode==None:
            download = ['download_chandra_obsid', 
                        f'{self.obsids}', 
                        'asol,bias,bpix,eph1,evt1,flt,msk,mtl,pbk,stat']
            run(download)

            for obs in self.obsid_list:
                os.mkdir(f'./{obs}/reduction/')
                logging.info(f"Created target directory:\n{os.path.abspath(f'./{obs}/reduction/')}\n")
                os.mkdir(f'./{obs}/repro/')
                logging.info(f"Created target directory:\n{os.path.abspath(f'./{obs}/repro/')}\n")
            
            gz_files = []
            for root, dir, files in os.walk('.', topdown=True):
                for name in files:
                    if name.endswith('.gz'):
                        gz_files.append(os.path.join(root, name))
            if gz_files:
                logging.info('Gunzipping files...')
                for file_path in gz_files:
                    with gzip.open(file_path, 'rb') as f_in:
                        with open(file_path.replace('.gz', ''), 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(file_path)
                    logging.info(f'Gunzipped {file_path}.')
            else:
                logging.info('No files to gunzip!')
        elif mode=='q':
            pass
        
    def bary_corr(self, infile):
        """
        Checks orbit file and selects appropriate one if there are multiple. Applies
        barycenter corrections to given a file.
        ----------
        infile : obj
            File to apply barycenter corrections
        """
        orbits = [i for i in os.listdir('primary/') if i.endswith('eph1.fits')]
        if len(orbits) == 1:
            self.orbit_file = os.path.abspath(f'primary/{orbits[0]}')
        elif len(orbits) > 1:
            times = []
            for i in orbits:
                if i.startswith('orbitf0'):
                    times.append(i.split('orbitf0')[1].split('N')[0])
                else:
                    times.append(i.split('orbitf')[1].split('N')[0])
                    
            tstart = fits.open(self.evt1_file)[1].header['TSTART']
            valid_times = [i for i in times if float(i) < tstart]
            chosen = max(valid_times)
            orbit_file = [i for i in orbits if chosen in i][0]
            self.orbit_file = os.path.abspath(f'primary/{orbit_file}')

        axbary.punlearn()
        logging.info(f'Applying barycenter correction to {infile} file...')
        if infile == self.asol_file:
            self.bary_asol = infile.replace('asol1', 'bary_asol1').replace('reduction', 'repro')
            axbary(infile=infile,
                   orbitfile=self.orbit_file,
                   outfile=self.bary_asol,
                   ra=self.RA, dec=self.DEC,
                   clobber='yes')

        elif infile == self.repro_evt2:
            self.bary_evt2 = infile.replace('repro', 'bary').replace('reduction', 'repro')
            axbary(infile=infile,
                   orbitfile=self.orbit_file,
                   outfile=self.bary_evt2,
                   ra=self.RA, dec=self.DEC,
                   clobber='yes')
            
        elif infile == self.stat_file:
            self.bary_stat = infile.replace('stat1.fits', 'bary_stat1.expstats').replace('reduction', 'repro')
            axbary(infile=infile,
                   orbitfile=self.orbit_file,
                   outfile=self.bary_stat,
                   ra=self.RA, dec=self.DEC,
                   clobber='yes')

            
    def reprocess(self):
        """
        Reprocess data for each observation, bary corrects the asol, reprocessed
        event, & stat files, and updates the corresponding event file headers 
        """
        for obs in self.obsid_list:
            os.chdir(f'{obs}/')
            event = [i for i in os.listdir('secondary/') if 'evt1' in i][0]
            self.evt1_file = os.path.abspath(f'secondary/{event}')
            
            self.read_mode()

            logging.info(f'Reprocessing data for obsid {obs}...')
            ardlib.punlearn()
            repro_func = ['chandra_repro',
                          f'indir=./',
                          f'outdir=reduction/',
                          'clobber=yes',
                          'verbose=1']
            self.verbose(repro_func)

            os.chdir('reduction/')
            repro_files = ['repro_bpix1.fits', 'repro_flt2.fits', 'repro_fov1.fits', 'msk1.fits']

            for name in os.listdir():
                if name.endswith('repro_evt2.fits'):
                    self.repro_evt2 = os.path.abspath(name)
                elif name.endswith('stat1.fits'):
                    self.stat_file = os.path.abspath(name)
                elif name.endswith('asol1.fits'):
                    self.asol_file = os.path.abspath(name)
                for file in repro_files:
                    if name.endswith(file):
                        shutil.copy2(name, '../repro/')
                        path = os.path.abspath(name).replace('reduction', 'repro')
                        setattr(self, file[:-6], path)

            os.chdir('../')
            self.bary_corr(self.asol_file)
            self.bary_corr(self.stat_file)
            self.bary_corr(self.repro_evt2)
            
            keys = ['ASOLFILE', 'BPIXFILE', 'FLTFILE', 'MASKFILE']
            values = [self.bary_asol, self.repro_bpix, self.repro_flt, self.msk]
            for key, value in zip(keys, values):
                self.edit_file(self.bary_evt2, key, value)
            shutil.rmtree('primary/')
            shutil.rmtree('secondary/')
            
            os.chdir('../')


    def merge(self):
        """
        Merges all reprocessed observations and view target in DS9
        """
        if os.path.exists('multi_obs/') == True:
            pass
        elif os.path.exists('multi_obs/') == False:
            os.mkdir('multi_obs/')
            logging.info(f"Created target directory:\n{os.path.abspath('multi_obs/')}\n")

            for root, dirs, files in os.walk('./', topdown=True):
                for name in files:
                    if name.endswith('bary_evt2.fits'):
                        if 'multi_obs' not in root:
                            path = os.path.join(root, name)
                            shutil.copy2(path, 'multi_obs/')

            if self.continuous:
                for obs in self.continuous:
                    self.obsid_list.remove(obs)
        
            logging.info('Merging observations...')
            merge_obs.punlearn()
            merge = ['merge_obs',
                     f'infiles={self.obsids}',
                     'outroot=multi_obs/',
                     'binsize=2',
                     f'bands={self.total}',
                     'clobber=yes',
                     'verbose=1']
            self.verbose(merge)

        os.chdir('multi_obs/')


    def periodogram(self, fmax=1e-2, step=1e-6, figure=None):
        """
        Produces a periodogram using the Rayleigh statistic
        and identifies the strongest signal
        ----------
        times : array
            Target event times
        fmin, fmax, step : float, optional
            The desired frequency range and step size
        figure : None, 'y'
            Switch to save figure
        """
        times = read_file('merged_evt.fits[sky=region(source.reg)]').get_column('TIME').values

        gti = []
        obs_file = [i for i in os.listdir() if 'bary' in i]
        for i in obs_file:
            hdul = fits.open(i)
            hdr = hdul[1].header
            gti.append(hdr['TSTOP'] - hdr['TSTART'])

        if self.continuous:
            continuous = []
            for obs in self.continuous:
                for i in os.listdir():
                    if i in obs:
                        continuous.append(i)
            for c, obs in zip(continuous, self.continuous):
                self.position(c)
                self.ds9_viewer(c, self.RA, self.DEC, reg=f'region_{obs}', img=f'ds9_{obs}', src=f'source_{obs}', bkg=f'background_{obs}')
                ctimes = read_file(f'{c}[sky=region(source_{obs}.reg)]').get_column('TIME').values
                times.append(ctimes)

        N = len(times)
        fmin = 1/min(gti)
        if fmin >= 5e-4:
            fmin = 1e-5

        frequencies = np.arange(fmin, fmax, step)
        fundamental = np.zeros_like(frequencies)
        harmonic = np.zeros_like(frequencies)

        for idx, freq in enumerate(frequencies):
            phi = 2 * np.pi * freq * times
            fundamental[idx] = (2 / N) * (np.sum(np.cos(phi))**2 + np.sum(np.sin(phi))**2)
            harmonic[idx] = fundamental[idx] + (2 / N) * (np.sum(np.cos(2 * phi))**2 + np.sum(np.sin(2 * phi))**2)

        fig, ax = plt.subplots(2, 1, figsize=(12, 4), sharex=True)
        fig.subplots_adjust(hspace=0)

        for i, pow in enumerate([fundamental, harmonic]):
            peaks, _ = find_peaks(pow)
            peak_idx = peaks[np.argmax(pow[peaks])]
            fpeak = frequencies[peak_idx]
            ffpeak = 1 / fpeak

            ax[i].plot(frequencies, pow, lw=0.9, color='black', label=f'$Z^2_{i+1}$')
            ax[i].axvline(fpeak, ls='dashed', lw=0.8, color='red', label=f'{ffpeak:.2f} s')
            ax[i].set_ylim(bottom=0)
            ax[i].legend(fontsize=10.75)
        
        ax[0].set_xlim(fmin, fmax)
        labels = [tick.get_text() for tick in ax[0].get_yticklabels()]
        labels[0] = ''
        ax[0].set_yticks(ax[0].get_yticks())
        ax[0].set_yticklabels(labels, fontsize=11)

        ax[1].tick_params(axis='both', labelsize=11)
        ax[1].set_xlabel('Frequency [Hz]', fontsize=12, labelpad=10)

        fig.supylabel('Power', x=0.07, fontsize=13)
        fig.suptitle(f'Rayleigh Statistic Power Spectrum of {self.target}', y=0.95, fontsize=14)

        if figure == None:
            plt.show()
        elif figure == 'y':
            plt.savefig('power_spectrum.eps')   
            plt.show()

