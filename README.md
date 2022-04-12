# Dust maps in the Magellanic Clouds

We present high-resolution maps of the dust reddening in the Magellanic Clouds (MCs). The maps cover the Large and Small Magellanic Cloud (LMC and SMC) area and have a spatial angular resolution between ~ 26 arcsec and 55 arcmin. Based on the data from the optical and near-infrared (IR) photometric surveys, including the Gaia Survey, the SkyMapper Southern Survey (SMSS), the Survey of the Magellanic Stellar History (SMASH), the Two Micron All Sky Survey (2MASS) and the near-infrared YJKS VISTA survey of the Magellanic Clouds system (VMC), we have obtained multi-band photometric stellar samples containing over 6 million stars in the LMC and SMC area. Based on the measurements of the proper motions and parallaxes of the individual stars from Gaia Early Data Release 3 (Gaia EDR3), we have built clean samples that contain stars from the LMC, SMC and Milky Way (MW), respectively. We apply the spectral energy distribution (SED) fitting to the individual sample stars to estimate their reddening values. As a result, we have derived the best-fitting reddening values of ~ 1.9 million stars in the LMC, 1.5 million stars in the SMC and 0.6 million stars in the MW, which are used to construct dust reddening maps in the MCs. Our maps are consistent with those from the literature. The resultant high-resolution dust maps in the MCs are not only important tools for reddening.

For more details, one can refer to our published paper at https://ui.adsabs.harvard.edu/abs/2022MNRAS.tmp..102C/abstract.

# Data access

Our reddening maps can be accessed at: http://paperdata.china-vo.org/diskec/mcdust/mcmap.fits.

The final catalogue containing the best-fit values of E(B-V) from our SED fitting of about 4 million stars can be accessed at: http://paperdata.china-vo.org/diskec/mcdust/starebv.fits. 

The Milky Way foreground reddening maps of the area toward the Magellanic Clouds can be accessed at: http://paperdata.china-vo.org/diskec/mcdust/mwmap.fits. 

# How to use the maps

To use our maps, we have provided a pyton procedure 'chen2021.py'. The procedure relies on the `dustmaps` package (https://github.com/gregreen/dustmaps).

The following steps show how to install our procedure:
1. Install the `dustmaps` package (pip install dustmaps).
2. Download the procedure 'chen2021.py' in this project and put the procedure in the `dustmaps` directory. 
3. Download the reddening maps (http://paperdata.china-vo.org/diskec/mcdust/mcmap.fits) to the '/path/where/you/want/large/data/files/stored'.


An example is given below to show how to obtain E(B-V) values from our dust maps:

    >>> from dustmaps.config import config
    >>> config['data_dir'] = '/path/where/you/want/large/data/files/stored'
    >>> from dustmaps.chen2021 import Chen2021Query
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>>
    >>> ebv = Chen2021Query()
    >>> ra = np.array([81.28,12.80])
    >>> dec = np.array([-69.78,-73.15])
    >>> c = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')  
    >>> ebv(c)
        array([[0.1063],[0.0076]], dtype=float32) 
     
--------------------------------------------------------------------------------
If you have any questions or need more informations, please send emall to Bingqiu Chen (bchen@ynu.edu.cn) and Helong Guo (helong_guo@mail.ynu.edu.cn).

--------------------------------------------------------------------------------

# Conditions for using the program

The program relies on the [`astropy.coordinates.SkyCoord`](http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord) and [`dustmaps`](https://github.com/gregreen/dustmaps)package. 
The `dustmaps` package provides a uniform interface for dealing with a number of 2D and 3D maps of interstellar dust reddening/extinction. For details and how to install, please refer to https://github.com/gregreen/dustmaps
