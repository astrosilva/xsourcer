import logging
from subprocess import Popen, PIPE, run
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
from astropy.io import fits, ascii
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
        self.target = target.replace(' ', '_')
        self.name = self.target.split('_')[1]
        self.ra = None
        self.dec = None
        self.RA = None
        self.DEC = None
        self.obsids = None
        self.obsid_list = None
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
        infile : obj
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
        infile : obj
            The original file to by copied
        outfile : str, obj
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


    def position(self):
        """
        Read the target RA and Dec for a given observation
        """
        hdul = fits.open(self.evt1_file)
        hdr = hdul[1].header
        self.RA, self.DEC = str(hdr['RA_TARG']), str(hdr['DEC_TARG'])
        logging.info(f'RA: {self.RA},  DEC: {self.DEC}')

        
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
            self.position()
            thresh = 0.5 / 3600 
            RA_frac = abs(self.RA - round(self.RA, 6))
            DEC_frac = abs(self.DEC - round(self.DEC, 6))

            if RA_frac < thresh or DEC_frac < thresh:
                logging.error(f'Insufficient coordinate precision! Cross-reference with another catalog and manually edit RA_TARG and DEC_TARG in {self.evt1_file}')
                exit()
            else:
                logging.info('Adequate CC coordinate precision.')
    
    def ds9_viewer(self, infile):
        logging.info('Loading target...')
        ds9_view = ['ds9',
                    infile,
                    '-scale', 'log',
                    '-colorbar', 'no',
                    '-cmap', '2', '0.4',
                    '-bin', 'buffersize', '8192',
                    '-bin', 'factor', '2',
                    '-pan', 'to', self.ra, self.dec, 'wcs', 'icrs',
                    '-zoom', '3',
                    '-region', 'format', 'ds9',
                    '-region', 'system', 'wcs',
                    '-region', 'sky', 'icrs',
                    '-region', 'skyformat', 'sexagesimal',
                    '-region', 'command', f'circle {self.ra} {self.dec} 4',
                    '-region', 'command', f'annulus {self.ra} {self.dec} 6 11',
                    '-region', 'system', 'physical',
                    '-region', 'format', 'ciao',
                    '-region', 'save', 'regions.reg',
                    '-saveimage', 'ds9_target.eps']
        run(ds9_view)
        input('Close DS9 window and press enter when finished. ')
        

        regions = ascii.read('regions.reg', format='no_header')
        ascii.write(regions[0], 'source.reg', format='no_header', overwrite=True)
        ascii.write(regions[1], 'background.reg', format='no_header', overwrite=True)