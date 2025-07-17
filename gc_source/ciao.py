# Fundamentals
import os
import gzip
import shutil
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
from subprocess import run

# Astropacks
from astroquery.simbad import Simbad
from PyAstronomy.pyasl import coordsDegToSexa as DtS
from astropy.io import ascii, fits
from pycrates import read_file
from ciao_contrib.runtool import *
from func import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

obsid = '7036'
os.chdir(obsid)

initialize = True
reprocess = True
bary_correct = True
source_extract = False

groupsize, binsize = 2,	5
group_spectrum = False
lightcurve = False


if initialize == True:
    hidden = []
    for root, directory, files in os.walk('.', topdown=True):
        for name in files:
            if name.startswith('._'):
                hidden.append(os.path.join(root, name))
    if hidden:
        logging.info('Removing hidden files...')
        for file_path in hidden:
            os.remove(file_path)
            logging.info(f'Deleted {file_path}')
    else:
        logging.info('No hidden files found.')
            

    gz_files = []
    for root, directory, files in os.walk('.', topdown=True):
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
            logging.info(f'Gunzipped {file_path}')
    else:
        logging.info('No files to gunzip.')
        
    
    if os.path.exists('reduction/') == True:
        logging.info('Reduction directory exists.')
        pass
    elif os.path.exists('reduction/') == False:
        logging.info('Creating reduction directory...')
        os.mkdir('reduction/')
    
    primary_links = ['evt2.fits', 'bpix1.fits', 'fov1.fits', 'eph1.fits', 'asol1.fits']
    os.chdir('primary/')
    primary_list = os.listdir()
    link_files = [os.path.join(os.getcwd(), i) for i in primary_list for j in primary_links if i.endswith(j)]

    secondary_links = ['evt1.fits', 'flt1.fits', 'msk1.fits', 'mtl1.fits', 'stat1.fits', 'pbk0.fits']
    os.chdir('../secondary/')
    secondary_list = os.listdir()
    link_files += [os.path.join(os.getcwd(), i) for i in secondary_list for j in secondary_links if i.endswith(j)]

    link_files = sorted(link_files)
    links = primary_links + secondary_links
    
    os.chdir('../reduction/')
    logging.info('Creating symbolic links...')
    for files, symbols in zip(link_files, links):
        dest = os.path.join(os.getcwd(), symbols)
        check_lf = files.split('_')[-1]
        
        if check_lf == symbols:
            try:
                os.symlink(files, dest)
            except:
                logging.info(f"{dest.rsplit('reduction/')[1]} already exists! Skipping...")
        elif check_lf != symbols:
            logging.error(f"ERROR: {dest.rsplit('reduction/')[1]} is not assigned to the correct file!")
            exit()

elif initialize == False:
    os.chdir('reduction')

    
if reprocess == True:
    readmode = ['dmkeypar',
                'infile=evt1.fits',
                'keyword=readmode',
                'echo=yes']
    mode = ''
    with Popen(readmode, stdout=PIPE, text=True) as mod:
        for line in mod.stdout:
            mode += line.replace('\n', '')

    if mode == 'CONTINUOUS':
        hdul = fits.open('evt1.fits')
        hdr = hdul[1].header
        RA, DEC = str(hdr['RA_TARG']), str(hdr['DEC_TARG'])
        print(RA)
        RA_dec = RA[RA.index('.') + 1:]
        DEC_dec = DEC[DEC.index('.') + 1:]

        if len(RA_dec) < 4 or len(DEC_dec) < 4:
            logging.error('Inadequate coordinate accuracy precision!')
            exit()
            #dmhedit acis_evt1.fits filelist="" operation=add key=RA_TARG value='newvalue'
            #dmhedit acis_evt1.fits filelist="" operation=add key=DEC_TARG value='newvalue'
        else:
            logging.info('Adequate coordinate accuracy precision.')

    logging.info('Reprocessing data...')
    ardlib.punlearn()
    repro_func = ['chandra_repro',
                  os.getcwd().replace('reduction', ''),
                  os.getcwd(),
                  'clobber=yes',
                  'verbose=1']
    verbose(repro_func)

elif reprocess == False:
    logging.info('Data reprocessing skipped.')
    
    
infile = 'spec_source.pi'
key_val = [('BACKFILE', infile.replace('.pi', '_bkg.pi')),
           ('RESPFILE', infile.replace('.pi', '.rmf')),
           ('ANCRFILE', infile.replace('.pi', '.corr.arf'))]

evt2 = str([i for i in os.listdir() if i.endswith('repro_evt2.fits')][0])
asol1 = str([i for i in os.listdir() if i.endswith('asol1.fits')][0])
evt2_bary = evt2.replace('repro', 'bary')
asol1_bary = asol1.replace('asol1', 'bary_asol1')


if bary_correct == True:
    logging.info('Applying barycenter correction...')
    
    hdul = fits.open(evt2)
    hdr = hdul[1].header
    RA = ''
    DEC = ''

    for key in hdr:
        if 'RA_TARG' in key:
            RA += str(hdr[key])
            print(f'{key} = {hdr[key]}')
        elif 'DEC_TARG' in key:
            DEC += str(hdr[key])
            print(f'{key} = {hdr[key]}')

    bary_corr(evt2, RA, DEC)
    bary_corr(asol1, RA, DEC)
    bary_corr('stat1.fits', RA, DEC)
    
    edit_file(evt2_bary, 'ASOLFILE', asol1_bary)

elif bary_correct == False:
    logging.info('Barycenter correction skipped.')


    
if source_extract == True:
    logging.info('Extracting source spectrum...')
    target = os.getcwd().split('/')[-3].replace('_', ' ')

    query = Simbad.query_object(target)
    ra = query['ra'].value
    dec = query['dec'].value
    sexa = DtS(ra, dec, fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f'))
    ra, dec = sexa.split(' ')

    ds9_view = ['ds9',
                evt2,
                '-scale', 'log',
                '-colorbar', 'no',
                '-cmap', '4', '0.25',
                '-bin', 'about', '8192', '8192',
                '-bin', 'about', 'center',
                '-bin', 'factor', '4',
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
                '-region', 'save', 'regions.reg']
    run(ds9_view)

    regions = ascii.read('regions.reg', format='no_header')
    source = ascii.write(regions[0], 'source.reg', format='no_header', overwrite=True)
    background = ascii.write(regions[1], 'bkg.reg', format='no_header', overwrite=True)
        
    specextract.punlearn()
    spec_extract = ['specextract',
                    f'infile={evt2_bary}[sky=region(source.reg)]',
                    f'bkgfile={evt2_bary}[sky=region(bkg.reg)]',
                    'outroot=spec_source',
                    'correctpsf=yes',
                    'weight=no',
                    'grouptype=NONE',
                    'binspec=NONE',
                    f'asp={asol1_bary}',
                    'dafile=CALDB',
                    'mskfile=msk1.fits',
                    'badpixfile=bpix1.fits',
                    'clobber=yes',
                    'verbose=1']
    verbose(spec_extract)

    for key, value in key_val:
        edit_file(infile, key, value)

elif source_extract == False:
    logging.info('Source extraction skipped.')

    
        
if group_spectrum == True:
    grouped = group_spec(infile, groupsize=groupsize)
    for key, value in key_val:
        edit_file(grouped, key, value)

        

if lightcurve == True:
    logging.info('Producing lightcurve...')
    evt2_src = read_file(f'{evt2_bary}[sky=region(source.reg)]')
    evt2_bkg = read_file(f'{evt2_bary}[sky=region(bkg.reg)]')
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
            
    light_curve(evt2_bary, ccd_val, binsize=binsize)



