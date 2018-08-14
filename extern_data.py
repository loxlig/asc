import ephem
from datetime import datetime
import numpy as np
from astropy.io import fits
from glob import glob


class Data(object):
    def __init__(self):
        self.obsname = 'LO'
        self.location = 'DCT'
        self.temperature = 19.3  # not known yet
        self.rel_humidity = 61  # not known yet
        self.exposure = 0
        self.date = ''
        self.time = ''
        self.lst = ''
        self.moon_illum = 0
        self.moon_elev = 0
        self.sun_elev = 0
        self.seqnum = 0
        self.heatstat = 'OFF'  # not known yet

    def external_data(self):
        lastfile = sorted(glob('/Users/Desktop/asc/TARGET__*.fit'))[-1]
        curfits = fits.open(lastfile)
        expt = curfits[0].header['EXPTIME']

        seqnum = lastfile[29:32]
        dct = ephem.Observer()
        dct.lat, dct.lon = '34.7443', '-111.4223'
        dct.elevation = 2361
        utcdate = curfits[0].header['DATE-OBS'][0:10]
        utctime = curfits[0].header['DATE-OBS'][11:19]
        dct.date = utcdate + ' ' + utctime
        lst = dct.sidereal_time()

        moon = ephem.Moon()
        moon.compute(dct.date, epoch=dct.date)
        moon_loc = ephem.Moon()
        moon_loc.compute(dct)
        moon_alt = (moon_loc.alt)
        moondms = np.degrees(ephem.degrees(moon_alt))
        moonflt = '{0:02.1f}'.format(moondms)
        mphase = (moon.phase) 
        mphaseflt = '{0:04.1f}'.format(mphase)
    
        sun = ephem.Sun()
        sun.compute(dct)
        sun_alt = (sun.alt)
        sundms = np.degrees(ephem.degrees(sun_alt))
        sunflt = '{0:02.1f}'.format(sundms)
    
        self.exposure = expt
        self.date = utcdate          
        self.time = utctime         
        self.lst = str(lst)[0:8]      
        self.moon_illum = mphaseflt
        self.moon_elev = moonflt
        self.sun_elev = sunflt
        self.seqnum = seqnum

        return self.obsname, self.location, self.temperature, self.rel_humidity, self.exposure, self.date, self.time, \
            self.lst, self.moon_illum, self.moon_elev, self.sun_elev, self.seqnum, self.heatstat


data = Data()
