# Fundamentals
import os
import gzip
import shutil
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
# Astropacks
from astroquery.simbad import Simbad
from PyAstronomy.pyasl import coordsDegToSexa as DtS
from astropy.io import ascii, fits
from pycrates import read_file
from ciao_contrib.runtool import *
from subprocess import run
from func import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

obsid = '2943'
os.chdir(obsid)

multi_corr = False
gaia_corr = False


if gaia_corr == True:
    logging.info('Applying astrometric corrections via Gaia...')
    dmcopy(infile='fov1.fits',
           outfile='fov1.fov',
           clobber='yes')
    
    dmcopy(infile='evt2.fits[sky=region(fov1.fov)]',
           outfile=f'evt2_filt.fits',
           clobber='yes')

    flux_img = ['fluximage',
                'infile=evt2_filt.fits',
                'outroot=./',
                'bands=broad',
                'binsize=1',
                'psfecf=0.5',
                'clobber=yes',
                'verbose=1']
    verbose(flux_img)

    wavdetect.punlearn()
    wavedetect = ['wavdetect',
                  'infile=broad_thresh.img',
                  'outfile=broad_src.fits',
                  'scell=.',
                  'imagefile=.',
                  'defnbkgfile=.',
                  'scales=1 2 4 6 8 12 16 24 32',
                  'expfile=broad_thresh.expmap',
                  'psffile=broad_thresh.psfmap',
                  'clobber=yes', 'verbose=1']
    verbose(wavedetect)

    ds9_view = ['ds9',
                'broad_thresh.img',
                '-scale', 'log',
                '-colorbar', 'no',
                '-cmap', '4', '0.5',
                '-zoom', 'to', 'fit',
                '-catalog', 'cds', 'I/355',
                '-catalog', 'symbol', 'color', 'red',
                '-region', 'broad_src.fits']
    run(ds9_view)

    wcs_match(infile='broad_src.fits',
              refsrcfile='ds9.tsv[opt skip=1][cols ra=col1,dec=col2]',
              outfile='crossmatch.xform',
              wcsfile='broad_thresh.img',
              method='trans',
              radius='5',
              clobber='yes')
    
    wcs_update(infile='asol1.fits',
               outfile='asol1_corr.fits',
               transformfile='crossmatch.xform',
               wcsfile='broad_thresh.img',
               clobber='yes')

    dmcopy(infile='evt2.fits',
           outfile='evt2_corr.fits',
           clobber='yes')
    
    wcs_update(infile='evt2_corr.fits',
               outfile='')

    edit_file('evt2_corr.fits', 'ASOLFILE', 'asol1_corr.fits')
    
elif gaia_corr == False:
    logging.info('Gaia astrometric correction skipped.')
    

