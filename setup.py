from distutils.core import setup

setup(
    name = "GammaLike",
    packages = ["GammaLike",], 
    #install_requires = ["numpy","scipy","h5py","pyfits","iminuit","healpy","cPickle","astropy"],
    version = "1.0.0",
    description = "Gmma-ray Template Regression Toolkit",
    author = "Eric Carlson",
    author_email = "eric.carlson.128@gmail.com",
    url = "http://planck.ucsc.edu/GammaLike/",
    download_url = "http://chardet.feedparser.org/download/python3-chardet-1.0.1.tgz",
    keywords = ["Gamma ray","dark matter","Fermi"],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Data Analysis",
        ]
)
