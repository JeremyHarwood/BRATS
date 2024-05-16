# Broadband Radio Astronomy Tools (BRATS)

BRATS is a software package for Linux and Mac which provides a variety of tools for the spectral analysis of both new broad-bandwidth radio data, along with legacy support for narrowband telescopes. Key features include:

- Adaptive Regions - Automatic selection of regions based on user parameters (e.g. signal to noise).
- Spectral Age Model Fitting - Fit models of spectral ageing on small spatial scales.
- Injection Index Minimisation - Automatic determination of the best-fitting injection index.
- Statistical Testing - Chi-squared testing, error maps, confidence levels and binning of model fits.
- Spectral Index Mapping - Map spectral index as a function of position.
- Image Reconstruction - Reconstruct sources at any frequency for a given model and parameter set.
- Image Subtraction - Subtract any two FITS images and output residual maps.
- Multi-frequency Image Combination - Easily combine and scale FITS images in the image plane.
- Resize radio maps - Quickly and easily crop or expand FITS images.

For a full list of features, usage instructions, and a walkthrough tutorial, please refer to the BRATS cookbook available on the GitHub repository.

**If you have made use of this software please cite [Harwood et al., 2013, MNRAS, 435, 3353](http://mnras.oxfordjournals.org/content/435/4/3353 "Spectral ageing in the lobes of FR-II radio galaxies: New methods of analysis for broadband radio data") and [Harwood et al., 2015, MNRAS, 454, 3403](http://mnras.oxfordjournals.org/content/454/4/3403 "Spectral ageing in the lobes of cluster-centre FR-II radio galaxies"). If possible, also include a link to this GitHub respository.**

The original papers detailing BRATS can be obtained on [MNRAS](http://mnras.oxfordjournals.org/content/435/4/3353 "Spectral ageing in the lobes of FR-II radio galaxies: New methods of analysis for broadband radio data") or [arXiv](http://arxiv.org/abs/1308.4137 "Spectral ageing in the lobes of FR-II radio galaxies: New methods of analysis for broadband radio data"), with more recent features detailed [here (MNRAS)](http://mnras.oxfordjournals.org/content/454/4/3403 "Spectral ageing in the lobes of cluster-centre FR-II radio galaxies") and [here (arXiv)](http://arxiv.org/abs/1509.06757v1 "Spectral ageing in the lobes of cluster-centre FR-II radio galaxies"). Useful information and testing of some features against simulated data can also be found [here (arXiv)](http://arxiv.org/abs/1409.1579v1 "Spectral age modelling of the Sausage cluster radio relic").

A Python wrapper and tools for automated image alignment can be found at: 

[BRATS Python Wrapper](https://github.com/JeremyHarwood/bratswrapper "BRATS Python Wrapper")

[Image Alignment Tool (DS9)](https://github.com/JeremyHarwood/bratsimagealignment "Image Alignment Tool (DS9)")

[Image Alignment Tool (CASA)](https://github.com/JeremyHarwood/bratsimagealigner_casa "Image Alignment Tool (CASA)")


## Known issues

- If the fits header is too long, the funtools library will cause the program to crash out with the error 'no WCS information in file while parsing filter at: XX:XX:XX.XXX'. This normally occurs due to a large history accumulated during reduction and can be fixed by running the AIPS Stalin command or history=False when exporting in CASA.
- There is a known bug in some versions of GSL that causes the chi-squared CDF function to fail, particularly for a large number of DoF or very low chi-squared values. The suppressconf command has been included as a work around but this can also sometimes throw up oddities when mapping later on. In such cases exporting the data via exportdata and mapping in an external program is recommended or, ideally, update GSL to a version where this issue has been resolved.
- When outputting in FITS format, certain header parameters (e.g. CTYPE) are grouped together rather than set out in the traditional format. This appears to be a FUNTOOLS issue.
- When using linear plots, the colour wedge scale is limited to 4 orders of magnitude. This appears to be a PGPLOT issue. The default is currently set to using a log scale which avoid this is most cases.
