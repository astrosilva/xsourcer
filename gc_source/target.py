# FUNCTIONALS
import logging
from subprocess import run
import numpy as np
import gzip, os, shutil
import matplotlib.pyplot as plt
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
# ASTROPACKS
from astropy.io import ascii, fits
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
        """
        if mode == None:
            if os.path.exists(f'{self.target}/') == True:
                os.chdir(f'{self.target}/')

            elif os.path.exists(f'{self.target}/') == False:
                os.mkdir(f'{self.target}/')
                logging.info(f'Created target directory:\n{os.getcwd()}\n')
                os.chdir(f'{self.target}/')

            self.root = os.getcwd()

            logging.info('Executing query...')
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

            ftable = ascii.read('query.dat')
            ftable = ftable[ftable['name'] == f' 2CXO {self.name}']
            ftable = ftable[ftable['instrument'] != ' HRC']
            ftable = ftable['name', 'obsid', 'ra', 'dec']
            ascii.write(ftable, 'query.dat', overwrite=True)
            print(ftable)

        elif mode == 'q':
            os.chdir(f'{self.target}/')
            ftable = ascii.read('query.dat')
        
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
            ex: download(254, 13615, 13176)
        retrieve : optional
            Switch for downloading data products. If False, only
            passes obsids
        """
        self.obsid_list = obsid_list
        self.obsids = ','.join(map(str, self.obsid_list))

        if mode == None:
            download = ['download_chandra_obsid', self.obsids]
            run(download)

            for obs in self.obsid_list:
                os.mkdir(f'./{obs}/reduction/')
                os.mkdir(f'./{obs}/repro/')
            
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
        
        elif mode == 'q':
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
            
            self.position()
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

            os.chdir('../')


    def merge_obs(self):
        """
        Merges all reprocessed observations and view target in DS9
        """
        if os.path.exists('multi_obs/') == True:
            pass
        elif os.path.exists('multi_obs/') == False:
            os.mkdir('multi_obs/')
            logging.info('Created multi_obs directory.')

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
        self.ds9_viewer('merged_evt.fits')
        
    
    def rayleigh(self, times, fmin=1e-5, fmax=1e-2, step=1e-6):
        N = len(times)
        frequencies = np.arange(fmin, fmax, step)
    
        fundamental = np.zeros_like(frequencies)
        harmonic = np.zeros_like(frequencies)
        for idx, freq in enumerate(frequencies):
            phi = 2 * np.pi * freq * times
            fundamental[idx] = (2 / N) * (np.sum(np.cos(phi))**2 + np.sum(np.sin(phi))**2)
            harmonic[idx] = fundamental[idx] + (2 / N) * (np.sum(np.cos(2 * phi))**2 + np.sum(np.sin(2 * phi))**2)

        return frequencies, fundamental, harmonic
    

    def periodogram(self):
        time = read_file('merged_evt.fits[sky=region(source.reg)]').get_column('TIME').values
        freq, fund, harm = self.rayleigh(time)

        fpeak =	freq[np.where(fund == max(fund))][0]
        ffpeak = 1 / fpeak

        hpeak = freq[np.where(harm == max(harm))][0]
        hhpeak = 1 / hpeak

        # Plot
        fig, ax = plt.subplots(2, 1, figsize=(12, 4), sharex=True)
        fig.subplots_adjust(hspace=0)

        ax[0].plot(freq, fund, lw=0.9, color='black', label=r'$Z^2_1$ (Fundamental)')
        ax[0].axvline(fpeak, ls='dashed', lw=0.8, color='red', label=f'{ffpeak:.2f} s')
        ax[1].plot(freq, harm, lw=0.9, color='black', label=r'$Z^2_2$ (Harmonic)')
        ax[1].axvline(hpeak, ls='dashed', lw=0.8, color='red', label=f'{hhpeak:.2f} s')
        
        ax[0].set_xlim(min(freq), max(freq))
        ax[0].set_ylim(bottom=0)
        labels = [tick.get_text() for tick in ax[0].get_yticklabels()]
        labels[0] = ''
        ax[0].set_yticks(ax[0].get_yticks())
        ax[0].set_yticklabels(labels)
        ax[0].legend()

        ax[1].set_ylim(bottom=0)
        ax[1].set_xlabel('Frequency [Hz]', fontsize=12, labelpad=5)
        ax[1].legend()

        fig.supylabel('Power', x=0.07, fontsize=12)
        fig.suptitle(f'Rayleigh Statistic Power Spectrum of {self.target}', y=0.95, fontsize=14)

        plt.savefig('power_spectrum.eps')
        plt.show()