xsourcer
=========

A tool for reprocessing and analyzing X-ray data for likely source identification.

---
Dependencies
--

Requires CIAO, HEASOFT recommended

https://cxc.cfa.harvard.edu/ciao/threads/ciao_install_conda
https://heasarc.gsfc.nasa.gov/docs/software/conda.html

---
`Using xsourcer`. 
--

Initialize the program and provide a target. 

Note: Currently, `xsourcer` is limited to Chandra data and only accepts targets in the Chandra-Discovered Source Naming Convention. 

```python
from target import Target
chandra = Target('CXOGC J174534.5-290201')
```

<br /> `query` : returns a list of all complete observations where the target is identified in the Chandra Source Catalog. <br />

Note: To mute `query` after the first run, pass `query(mode='q')`.
```python
chandra.query()
```

<br /> `download` : provide desired OBSIDs to download the required processing files. <br />

Note: To prevent redownloading data in subsequent runs, pass `mode='q'` after the last OBSID.
```python
chandra.download(2951, 3392, 13016, 22707)
```

<br /> `reprocess` : applies the corrections outlined by CIAO's ACIS reprocessing pipeline. <br />

Note: Once data is reprocessed, this task is no longer needed and can be commented out.
```python
chandra.reprocess()
```

<br /> Stack the reprocessed data with `merge`. Use `ds9_viewer` to verify the resulting image and extraction regions. <br />
```python
chandra.merge()
chandra.ds9_viewer('merged_evt.fits', chandra.ra, chandra.dec)
```

<br /> `periodogram` : generates a Rayleigh statistic power spectrum. <br />
```python
chandra.periodogram()
```
