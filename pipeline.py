# FUNCTIONALS
import numpy as np
import gzip, logging, os, shutil, sys
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from subprocess import Popen, PIPE, run
# ASTROPACKS
from astropy.io import fits
from astropy.table import unique, Table
from astroquery.simbad import Simbad
from pycrates import read_file
from ciao_contrib.runtool import *
from PyAstronomy.pyasl import coordsDegToSexa as DtS

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

"""
Parent class containing fundamental tasks and initial objects
"""

class Source:
    """
    target : str
        Should be of the format 2CXO/CXOGC JHHMMSS.s-DDMMSS
    """
    def __init__(self, target):
        self.target = target
        self.target_path = target.replace(' ', '_')
        self.name = target.split(' ')[1]
        self.ra = None
        self.dec = None
        self.RA = None
        self.DEC = None
        self.obsid_list = None
        self.continuous = []
        self.bands = [('broad', '0.5-7.0', '500:7000'),
                      ('soft', '0.5-2.0', '500:2000'),
                      ('hard', '2.0-7.0', '2000:7000')]
        
        self.asol_file = None
        self.evt1_file = None
        self.msk_file = None
        self.orbit_file = None
        
        self.repro_bpix = None
        self.repro_evt2 = None
        self.repro_flt = None
        self.repro_fov = None
        self.repro_stat = None
                
        self.bary_asol = None
        self.bary_evt2 = None
        self.bary_stat = None

        self.fov_fov = None
        self.filt_evt2 = None
        self.corr_asol = None
        self.corr_evt2 = None


    @property
    def obsids(self):
        return ','.join(map(str, self.obsid_list))
    
    
    def verbose(self, module):
        """
        Print status and progress for a given CIAO task
        ----------
        module : str
            A CIAO command line argument defined in a subsequent function
        """
        with Popen(module, stdout=PIPE, text=True) as task:
            for line in task.stdout:
                print(line, end='')

                
    def edit_file(self, infile, key, key_file):
        """
        Edit headers for a given FITS file
        ----------
        infile : str
            Typically the preferred events file
        key : str
            The desired header to edit
        key_file : str
            Path to desired update file
        """
        edit = ['dmhedit',
                f'infile={infile}',
                'filelist=NONE',
                'operation=add',
                f'key={key}',
                f"value='{key_file}'",
                'verbose=1']
        self.verbose(edit)


    def copy_file(self, infile, outfile, dir_name='.'):
        """
        Copy original file or copy and apply filter/binning
        ----------
        infile : str
            The original file to by copied
        outfile : str
            Name of new file
        dir_name : str, optional
            Path to desired directory
        """
        dmcopy.punlearn()
        copy = ['dmcopy',
                f'infile={infile}',
                f'outfile={dir_name}/{outfile}',
                'clobber=yes',
                'verbose=1']
        self.verbose(copy)


    def position(self, infile):
        """
        Read the target RA and Dec for a given observation
        """
        hdul = fits.open(infile)
        hdr = hdul[1].header
        self.RA, self.DEC = hdr['RA_TARG'], hdr['DEC_TARG']


    def read_mode(self):
        """
        Check if ACIS observation used Continuous Clocking mode and,
        if so, verify the coordinates are of the appropriate precision
        """
        readmode = ['dmkeypar',
                    f'infile={self.evt1_file}',
                    'keyword=readmode',
                    'echo=yes']
        mode = ''
        with Popen(readmode, stdout=PIPE, text=True) as task:
            for line in task.stdout:
                mode += line.replace('\n', '')
        
        if mode == 'CONTINUOUS':
            obsid = self.evt1_file[5:9]
            self.continuous.append(obsid)

            self.position(self.evt1_file)
            thresh = 4
            RA_frac = str(self.RA).split('.')[-1].rstrip('0')
            DEC_frac = str(self.DEC).split('.')[-1].rstrip('0')
            
            if len(RA_frac) < thresh or len(DEC_frac) < thresh:
                logging.error(f'Insufficient coordinate precision! Cross-referencing with SIMBAD\nRA: {self.RA}\nDec: {self.DEC}')
                search = Simbad.query_object(self.target)
                new_RA = search['ra'].value[0]
                new_DEC = search['dec'].value[0]

                if len(new_RA_frac) < thresh or len(new_DEC_frac) < thresh:
                    logging.error('Insufficient coordinate precision! Check another catalog.')
                    print(f'RA: {new_RA}\nDec: {new_DEC}')
                else:
                    logging.info('Updating RA and DEC...')
                    self.edit_file(self.evt1_file, 'RA_TARG', new_RA)
                    self.edit_file(self.evt1_file, 'DEC_TARG', new_DEC)
            else:
                logging.info('Adequate CC coordinate precision.')
        

    def ds9_viewer(self, infile, ra, dec, reg='regions', img='ds9_target', src='source', bkg='background'):
        """
        Opens FITS image of desired target in ds9 and generates region files based on query coordinates
        ----------
        infile : str
            The desired image file
        ra, dec : str
            RA and Dec of target
        reg : str
            Output name for region file
        img : str
            Output name for ds9 capture
        src : str
            Output name for source region file
        bkg : str
            Output name for background region file          
        """
        logging.info('Loading target...')

        ds9_view = ['ds9',
                    infile,
                    '-scale', 'log',
                    '-colorbar', 'no',
                    '-cmap', '2', '0.4',
                    '-bin', 'buffersize', '8192',
                    '-bin', 'factor', '2',
                    '-pan', 'to', ra, dec, 'wcs', 'icrs',
                    '-zoom', '3',
                    '-region', 'format', 'ds9',
                    '-region', 'system', 'wcs',
                    '-region', 'sky', 'icrs',
                    '-region', 'skyformat', 'sexagesimal',
                    '-region', 'command', f'circle {ra} {dec} 3',
#                    '-region', 'command', f'annulus {ra} {dec} 6 11',
                    '-region', 'format', 'ciao',
                    '-region', 'system', 'physical',
                    '-region', 'color', 'green',
                    '-region', 'load', 'merged_expmap_src.reg',
#                    '-region', 'save', f'{reg}.reg',
                    '-saveimage', f'{img}.eps']
#        run(ds9_view)

        print('Enter centroid sexagesmial coordinates')
#        self.ra = input('RA: ')
#        self.dec = input('Dec: ')
#        input('Close DS9 window and press enter when finished.')

#        regions = ascii.read(f'{reg}.reg', format='no_header')
#        ascii.write(regions[0], f'{src}.reg', format='no_header', overwrite=True)
#        ascii.write(regions[1], f'{bkg}.reg', format='no_header', overwrite=True)



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
        barycenter corrections to a given file.
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
                elif name.endswith('msk1.fits'):
                    self.msk_file = os.path.abspath(name)
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
            values = [self.bary_asol, self.repro_bpix, self.repro_flt, self.msk_file]
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
                    if name.endswith('bary_evt2.fits') or name.endswith('bary_asol1.fits') or name.endswith('msk1.fits') or name.endswith('repro_bpix1.fits'):
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
                     'bands=broad',
                     'binsize=1',
                     'psfecf=0.9',
                     'psfmerge=exptime',
                     'clobber=yes',
                     'verbose=1']
            self.verbose(merge)

        os.chdir('multi_obs/')

        logging.info('Building source catalogue with wavdetect...')
        wavdetect.punlearn()
        sources = ['wavdetect',
                   'infile=broad_thresh.img',
                   'psffile=broad_thresh.psfmap',
                   'expfile=broad_thresh.expmap',
                   'scales=2 4 8 16 24 32 48',
                   'outfile=merged_expmap_src.fits',
                   'scellfile=merged_expmap_scell.fits',
                   'imagefile=merged_expmap_imgfile.fits',
                   'defnbkgfile=merged_expmap_nbgd.fits',
                   'regfile=merged_expmap_src.reg',
                   'clobber=yes',
                   'verbose=1']
        self.verbose(sources)


    def spectrum(self):
        evt_list = []
        for root, dirs, files in os.walk('./', topdown=True):
            for name in files:
                    if name.endswith('bary_evt2.fits'):
                        evt_list.append(name)
    
        spectra = []
        for i, evt in enumerate(evt_list):
            srcflux.punlearn()
            sourceflux = ['srcflux',
                        f'infile={evt}',
                        f'pos={self.ra}, {self.dec}',
                        'outroot=./',
                        'bands=broad',
                        'regions=simple',
                        'psfmethod=marx',
                        'conf=0.68',
                        'clobber=yes',
                        'verbose=1']
            self.verbose(sourceflux)

            src = Table.read('0001_srcreg.fits', format='fits')
            src_x = str(src['X'][0])
            src_y = str(src['Y'][0])
            src_r = str(src['R'][0][0])
            with open('source.reg', 'w') as f:
                f.write(f'circle({src_x},{src_y},{src_r})')

            bkg = Table.read('0001_bkgreg.fits', format='fits')
            bkg_x = str(bkg['X'][0])
            bkg_y = str(bkg['Y'][0])
            bkg_rin = str(bkg['R'][1][0])
            bkg_rout = str(bkg['R'][0][0])
            with open('bkg.reg', 'w') as f:
                f.write(f'annulus({bkg_x},{bkg_y},{bkg_rin},{bkg_rout})')

            specextract.punlearn()
            spec_extract = ['specextract',
                            f'infile={evt}[sky=region(source.reg)]',
                            f'bkgfile={evt}[sky=region(bkg.reg)]',
                            f'outroot=spec_{i}',
                            'correctpsf=yes',
                            'weight=no',
                            'grouptype=NONE',
                            'binspec=NONE',
                            'dafile=CALDB',
                            'clobber=yes',
                            'verbose=1']
            self.verbose(spec_extract)
          
            key_val = [('BACKFILE', f'spec_{i}_bkg.pi'),
                       ('RESPFILE', f'spec_{i}.rmf'),
                       ('ANCRFILE', f'spec_{i}.corr.arf')]
            bkg_val = [('RESPFILE', f'spec_{i}_bkg.rmf'),
                       ('ANCRFILE', f'spec_{i}.corr.arf')]
                
            for key, value in bkg_val:
                self.edit_file(f'spec_{i}_bkg.pi', key, value)
            for key, value in key_val:
                self.edit_file(f'spec_{i}.pi', key, value)
            
            spectra.append(f'spec_{i}.pi')

        spectra_list = ",".join(map(str, spectra))
        combine_spectra.punlearn()
        combine = ['combine_spectra',
                   f'src_spectra={spectra_list}',
                   'outroot=combo',
                   'clobber=yes',
                   'verbose=1']
        self.verbose(combine)
  

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
