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
