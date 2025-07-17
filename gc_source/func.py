import os
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
from subprocess import Popen, PIPE
from ciao_contrib.runtool import *

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def verbose(module):
    with Popen(module, stdout=PIPE, text=True) as mod:
        for line in mod.stdout:
            print(line, end='')

            
def bary_corr(infile, RA, DEC, orbitfile='eph1.fits'):
    axbary.punlearn()
    logging.info(f'Applying barycenter correction to {infile} file...')
    if 'repro' in infile:
        axbary(infile=infile,
               orbitfile=orbitfile,
               outfile=infile.replace('repro', 'bary'),
               ra=RA, dec=DEC,
               clobber='yes')
    elif 'asol1' in infile:
        axbary(infile=infile,
               orbitfile=orbitfile,
               outfile=infile.replace('asol1', 'bary_asol1'),
               ra=RA, dec=DEC,
               clobber='yes')
            
    elif infile == 'stat1.fits' and infile != 'stat1.fits':
        axbary(infile=infile,
               orbitfile=orbitfile,
               outfile=infile.replace('stat1.fits', 'bary_stat1.expstats'),
               ra=RA, dec=DEC,
               clobber='yes')

        
def edit_file(infile, key, value):
    dmhedit.punlearn()    
    edit = ['dmhedit',
	    f'infile={infile}',
            'filelist=NONE',
            'operation=add',
            f'key={key}',
            f'value={value}',
            'verbose=1']
    verbose(edit)

            
def group_spec(infile, groupsize):
    logging.info('Grouping spectrum...')
    dmgroup.punlearn()
    outfile = f"{infile.replace('.pi', f'_gr{groupsize}.pi')}"
    group = ['dmgroup',
             f'infile={infile}',
             f'outfile={outfile}',
             'grouptype=NUM_CTS',
             f'grouptypeval={groupsize}',
             'binspec=',
             'xcolumn=channel',
             'ycolumn=counts',
             'clobber=yes',
             'verbose=1']
    verbose(group)

    return outfile


def light_curve(infile, ccd_val, binsize):
    bands = [('total', '0.1-10', '100:10000'),('soft', '0.1-2', '100:2000'), ('hard', '2-10', '2000:10000')]

    dmextract.punlearn()
    for band in bands:
        lcextract = ['dmextract',
                     f'infile={infile}[ccd_id={ccd_val}, sky=region(source.reg), energy={band[2]}][bin time=::{binsize}]',
                     f'outfile=lc_{band[0]}_{band[1]}.fits',
                     f'bkg={infile}[ccd_id={ccd_val}, sky=region(bkg.reg), energy={band[2]}]',
                     'opt=ltc1',
                     'clobber=yes',
                     'verbose=3']

        verbose(lcextract)

