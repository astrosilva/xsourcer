import logging
import sys
from subprocess import Popen, PIPE, run
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
from astropy.io import fits, ascii
from astroquery.simbad import Simbad
from ciao_contrib.runtool import *

"""
Parent class containing fundamental tasks and initial objects
"""

class Source:
    """
    target : str
        Should be of the format CXOGC JHHMMSS.s-DDMMSS
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
        self.total = '0.1:10:1.9'
        self.bands = [('total', '0.1-10', '100:10000'),
                      ('soft', '0.1-2', '100:2000'),
                      ('hard', '2-10', '2000:10000')]
        
        self.asol_file = None
        self.evt1_file = None
        self.msk = None
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
                    '-region', 'command', f'circle {ra} {dec} 4',
                    '-region', 'command', f'annulus {ra} {dec} 6 11',
                    '-region', 'system', 'physical',
                    '-region', 'format', 'ciao',
                    '-region', 'save', f'{reg}.reg',
                    '-saveimage', f'{img}.eps']
        run(ds9_view)
        input('Close DS9 window and press enter when finished. ')

        regions = ascii.read(f'{reg}.reg', format='no_header')
        ascii.write(regions[0], f'{src}.reg', format='no_header', overwrite=True)
        ascii.write(regions[1], f'{bkg}.reg', format='no_header', overwrite=True)