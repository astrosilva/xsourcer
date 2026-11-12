    def extract_spectrum(self, group=False, groupsize=5, outroot='pha_spec'):
        key_val = [('BACKFILE', f'{outroot}_bkg.pi'),
                   ('RESPFILE', f'{outroot}.rmf'),
                   ('ANCRFILE', f'{outroot}.corr.arf')]
        
        if group == False:
            specextract.punlearn()
            spec_extract = ['specextract',
                            f'infile={self.bary_evt2}[sky=region(source.reg)]',
                            f'bkgfile={self.bary_evt2}[sky=region(background.reg)]',
                            f'outroot={outroot}',
                            'correctpsf=yes',
                            'weight=no',
                            'grouptype=NONE',
                            'binspec=NONE',
                            f'asp={self.bary_asol}',
                            'dafile=CALDB',
                            'mskfile=msk1.fits',
                            'badpixfile=bpix1.fits',
                            'clobber=yes',
                            'verbose=1']
            self.verbose(spec_extract)

            for key, value in key_val:
                self.edit_file(f'{outroot}.pi', key, value)
            
        
        if group == True:
            outfile = f'{outroot}_gr{groupsize}.pi'
            
            logging.info('Grouping spectrum...')
            dmgroup.punlearn()
            group = ['dmgroup',
                     f'infile={outroot}.pi',
                     f'outfile={outfile}',
                     'grouptype=NUM_CTS',
                     f'grouptypeval={groupsize}',
                     'binspec=',
                     'xcolumn=channel',
                     'ycolumn=counts',
                     'clobber=yes',
                     'verbose=1']
            self.verbose(group)

            for key, value in key_val:
                self.edit_file(outfile, key, value)
        
