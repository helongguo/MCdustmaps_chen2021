# MCdustmaps_chen2021
Based on the data from the optical and near-infrared (IR) photometric surveys, including the Gaia Early Data Release 3 (Gaia eDR3), the SkyMapper Southern Survey DR2 (SMSS DR2), the Survey of the Magellanic Stellar History (SMASH), the Two Micron All Sky Survey (2MASS) and the near-infrared VISTA survey of the Magellanic Clouds system (VMC), we have obtained a multi-band photometric stellar sample of over 9 million stars in the LMC and SMC area. We applied the spectral energy distribution (SED) fitting to the individual sample stars to estimate their reddening values. As a result, we have derived the best-fitting reddening values of about 1.9 million stars in the LMC, 1.5 million stars in the SMC and 0.6 million stars in the MW, which are used to construct dust reddening maps in the MCs. The maps cover the Large and Small Magellanic Cloud (LMC and SMC) area and have a resolution of about 26 arcsec – 55 arcmin.

Our reddening maps can be accessed at: http://paperdata.china-vo.org/diskec/mcdust/mcmap.fits

The final catalogue containing the best-fit values of E(B-V) from our SED fitting of about 4 million stars can be accessed at: http://paperdata.china-vo.org/diskec/mcdust/starebv.fits

To get the color excess E(B-V) 

    >>> from dustmaps.config import config
    >>> config['data_dir'] = '/path/where/you/want/large/data/files/stored'
    >>> import dustmaps.chen2021
    >>> dustmaps.chen2021.fetch()
    >>> from dustmaps.chen2021 import Chen2021Query
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>>
    >>> ebv = Chen2021Query()
    >>>
    >>> c = SkyCoord(ra=80*u.degree,dec=-78.0*u.degree,frame='icrs')
        
    >>> ebv(c)
        0.1967
  
For details, please refer to https://github.com/gregreen/dustmaps  
    
--------------------------------------------------------------------------------
If you have any questions or need more informations, please send emall to Bingqiu Chen (bchen@ynu.edu.cn) and Helong Guo (helong_guo@mail.ynu.edu.cn).

--------------------------------------------------------------------------------

# Conditions for using the program
The program relies on the [`astropy.coordinates.SkyCoord`](http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord) and [`dustmaps`](https://github.com/gregreen/dustmaps)package. The `dustmaps` package provides a uniform interface for dealing with a number of 2D and 3D maps of interstellar dust reddening/extinction. For details and how to install, please refer to https://github.com/gregreen/dustmaps
