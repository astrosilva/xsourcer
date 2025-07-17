Descriptions
=
### astrometry.py
This module performs source matching and astrometric corrections using GAIA DR3 data. (Not yet integrated into the pipeline)

### ciao.py
This module is responsible for reprocessing ACIS products downloaded from the Chandra Data Archive. The following tasks can be performed:
initializing data environment, executing chandra_repro, applying barycenter corrections, and creating both spectral 
and light curve files. (Working on integrating with multi_obs.py and open access)

### func.py
This is the import module used by the rest of the programs in this repository. It houses the following CIAO tasks: verbose, barycenter 
corrections, file editing, grouping spectra, and producing light curves.

### light_curve.py
This module plots the extracted light curves, with time converted to JD. (Working on phase-folded light curves)

### multi_obs.py
This module can query the Chandra Source Catalog for a given source, concatenate the observations, and produce a composite light curve. (Working on integrating with ciao.py)

### period.py
This module generates periodograms/power spectra using Lomb-Scargle and the Rayleigh statistic $(Z^2_1)$. (Working on determining spin period)

---
