# GammaLike

The GammaLike package is devoted to easing multi-linear template regression on gamma-ray data using a binned likelihood approach.  While this library was designed to analyse Fermi-LAT data, it is easily extensible to any analysis using healpix maps.  Many tutorials are provided in the `tutorials` folder. 

The focus of this software is twofold.  First, efficiency.  The fermi science tools are prohibitively slow when working with large regions of interest.  GammaLike can perform full-sky likelihood fits to fermi data in only a few minutes at .25 degree resolution.  Second, is ease of use.  After supplying flux maps, GammaLike automatically performs convoutions against instrumental response functions and supports arbitrary binning, masking, and pixel weighting schemes.  Written in pure python, it is also easy to add new functionality.

Full documentation is available at xxxxxxxxxxxxxxxx


Features include:
- Fixed spectrum or bin-by-bin fitting in energy
- Built in support for a variety of gamma-ray structures such as the fermi bubbles or isotropic templates. 
- Fast point source generation tool
- Simple interfaces to galprop output files. 
- Multi-threading and GPU support in performance critical routines.
- Output models to HDF5 format. 
- Support for Arbitrary Masking and Pixel weights. 
- Can impose external uncertainties on components. 
- Generate dark matter templates including common DM profiles, translations, ellipticities. 

This library was developed originally to study both global fits to Fermi-LAT data, and to study the GeV Galactic center excess [arxiv:xxxx.xxxxxx]

