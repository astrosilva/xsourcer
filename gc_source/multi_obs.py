# Fundamentals
import os
import shutil
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
from subprocess import run
# Astropacks
from astroquery.simbad import Simbad
from PyAstronomy.pyasl import coordsDegToSexa as DtS
from astropy.io import ascii, fits
from astropy.table import Table
from ciao_contrib.runtool import *
from pycrates import read_file
from func import *

query = False
initialize = True
concatenate = True
inspect = True
lightcurve = False

if query == True:
#    search_csc(pos='CXOGC J174543.4-285841',
#               radius='0.5',
#               outfile='query.dat',
#               bands='broad',
#               filetype='evt',
#               catalog='latest',
#               clobber='yes')
#    
    table = ascii.read('query.dat')
    ftable = table[table['name'] == ' 2CXO J174543.4-285841']
    dtable = ftable[ftable['instrument'] != ' HRC']
    download = ['download_chandra_obsid', '7036']

    verbose(download)
    
if initialize == True:
    obs_list = ['2943', '3392_comp', '3665_comp', '7036', '242_comp']

    if os.path.exists('multi_obs/') == True:
        logging.info('multi_obs directory exists.')
        pass
    elif os.path.exists('multi_obs/') == False:
        logging.info('Creating multi_obs directory...')
        os.mkdir('multi_obs/')
        os.mkdir('multi_obs/repro/')
	
    for obs in obs_list:
        for root, directory, files in os.walk(f'./{obs}/reduction/'):
            for name in files:
                if 'repro' in name and 'evt2' not in name and 'asol1' not in name:
                    shutil.copy(os.path.join(root, name), './multi_obs/repro/')
                elif 'bary' in name:
                    shutil.copy(os.path.join(root, name), './multi_obs/repro/')
                elif name.endswith('msk1.fits') and name != 'msk1.fits':
                    shutil.copy(os.path.join(root, name), './multi_obs/repro/')

os.chdir('multi_obs/')

if concatenate == True:
    merge_obs.punlearn()
    merge = ['merge_obs',
             'infiles=repro/',
             'outroot=./',
             'asolfiles=',
             'binsize=2',
             'badpixfiles=',
             'maskfiles=',
             'bands=0.1:10:2',
             'clobber=yes',
             'verbose=1']

    verbose(merge)

evt2 = str([i for i in os.listdir() if i == 'merged_evt.fits'][0])
    
if inspect == True:
    target = os.getcwd().rsplit('/')[-2]
    query = Simbad.query_object(target)
    ra = query['ra'].value
    dec = query['dec'].value
    sexa = DtS(ra, dec, fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f'))
    ra, dec = sexa.split(' ')

    ds9_view = ['ds9',
                evt2,
                '-scale', 'log',
                '-colorbar', 'no',
                '-bin', 'about', '8192', '8192',
                '-bin', 'about', 'center',
                '-bin', 'factor', '1',
                '-pan', 'to', ra, dec, 'wcs', 'icrs',
                '-region', 'format', 'ds9',
                '-region', 'system', 'wcs',
                '-region', 'sky', 'icrs',
                '-region', 'skyformat', 'sexagesimal',
                '-region', 'command', f'circle {ra} {dec} 4',
                '-region', 'command', f'annulus {ra} {dec} 5 10',
                '-region', 'centroid',
                '-region', 'system', 'physical',
                '-region', 'format', 'ciao',
                '-region', 'save', 'regions.reg']

    run(ds9_view)
    
    regions = ascii.read('regions.reg', format='no_header')
    source = ascii.write(regions[0], 'source.reg', format='no_header', overwrite=True)
    background = ascii.write(regions[1], 'bkg.reg', format='no_header', overwrite=True)


if lightcurve == True:
    logging.info('Producing lightcurve...')
    evt2_src = read_file(f'{evt2}[sky=region(source.reg)]')
    evt2_bkg = read_file(f'{evt2}[sky=region(bkg.reg)]')
    ccd_src = evt2_src.get_column('ccd_id').values
    ccd_bkg = evt2_bkg.get_column('ccd_id').values

    ccd_val = ''

    if min(ccd_src) != min(ccd_bkg):
        logging.error('Values for ccd_id do not match!')
        print(f'Source | min: {min(ccd_src)}, max: {min(ccd_src)}, mean: {min(ccd_src)}')
        print(f'Background | min: {min(ccd_bkg)}, max: {min(ccd_bkg)}, mean: {min(ccd_bkg)}')
        exit()
    else:
        ccd_val += str(min(ccd_src))

    light_curve(evt2, ccd_val, binsize=5)
