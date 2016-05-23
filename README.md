# GammaLike

The GammaLike package is devoted to easing multi-linear template regression on gamma-ray data using a binned likelihood approach.  While this library was designed to analyse Fermi-LAT data, it is easily extensible to any analysis using healpix maps.  Many tutorials are provided in the `tutorials` folder. 

The focus of this software is twofold.  First, efficiency.  The fermi science tools are prohibitively slow when working with large regions of interest.  GammaLike can perform full-sky likelihood fits to fermi data in only a few minutes at .25 degree resolution.  Second, is ease of use.  After supplying flux maps, GammaLike automatically performs convoutions against instrumental response functions and supports arbitrary binning, masking, and pixel weighting schemes.  Written in pure python, it is also easy to add new functionality.

Package documentation is available at http://planck.ucsc.edu/gammalike/html/


Features include:
- Fixed spectrum or bin-by-bin fitting in energy
- Built in support for a variety of gamma-ray structures such as the fermi bubbles or isotropic templates. 
- Fast point source generation tool
- Simple interfaces to galprop output files. 
- Multi-threading support in performance critical regions.
- Output models to HDF5 format. 
- Support for Arbitrary Masking and Pixel weights. 
- Can impose external uncertainties on components. 
- Generate dark matter templates including common DM profiles, translations, ellipticities. 

This library was developed originally to study both global fits to Fermi-LAT data, and to study the GeV Galactic center excess [arxiv:1603.06584]

# Modified Galprop Codes 
The modified Galprop codes can be found at https://github.com/erccarls/galprop_diff_r2504

See the README of that repository for usage information and do not forget to cite the original Galprop code
http://galprop.sourceforge.net  and  http://galprop.stanford.edu

If you don't mind, cite us too: 
  http://arxiv.org/abs/1510.04698
  http://arxiv.org/abs/1603.06584

# New Galprop Gas Maps
To use our modified Galprop code you will need several additional files to be placed in the $galprop_home/FITS folder:  

The 3D gas maps for use in the Galprop code above can be found at:
https://www.dropbox.com/sh/lno27rn44ybu5pl/AACzUsBgG4AY1EClif2RGd8ya?dl=0 

The 9 galactocentric ring gas maps used for galprop can be found at:
https://www.dropbox.com/s/nu1iqvrj5uiu2al/galprop_rbands.tar.gz?dl=0

Finally, a conveinient script for generating Galdef files in python can be found at 
https://github.com/erccarls/hyades_scripts/blob/master/RunGalprop.py,
although one will want to modify this appropriately and remove the cluster submission code at the end. 


# New Diffuse Emission Models
The emission models from http://arxiv.org/abs/1603.06584 can be found at the following links:

Mod_A from Calore et al (but with the XCO profile fitted):
	https://www.dropbox.com/s/r48vjnmijbxx6xh/Mod_A_2D_XCO_P8.hdf5?dl=0

The f_H2=0.2 'Canonical' emission model with no wind: 
	https://www.dropbox.com/s/uee55f64klsmfjj/mod_s_46_XCO_P8_corrected.hdf5?dl=0

The best fitting f_H2=0.25 model with a 600 km/s radial wind at the Galactic center:
	


These files are stored in HDF5 format with the templates for each Galactic diffuse component stored at the path '/templates/'.  Further usage instructions can be found toward the bottom of the 'Fitting a Galprop Model.ipynb' tutorials.


If you need these in cartesian format, take a look at the script here:
	https://github.com/erccarls/GammaLike_dev/blob/master/Healpix2Mapcube_merged.py
and modify as needed. 




