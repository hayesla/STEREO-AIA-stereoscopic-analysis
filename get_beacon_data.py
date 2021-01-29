import sunpy.map 
from sunpy.util.scraper import Scraper
from sunpy.time import parse_time, TimeRange
from sunpy.net import Fido, attrs as a
from parfive import Downloader
import glob
from dateutil.relativedelta import relativedelta

def get_euvi_beacon_data(timerange, spacecraft="ahead", earliest=False):
    """
    Function to return url of EUVI beacon data for download.
    
    Parameters
    ----------
    timerange : ~`sunpy..time.TimeRange`
        timerange over which to do the search.
    spacecraft : ~`str`, default = `ahead`
        which STEREO spacecraft to use - either ahead or behind.
        
    Returns
    -------
    fileslist
    
    """
    spacecraft_letter = spacecraft.capitalize()[0]

    # define pattern
    stereo_pattern = ("https://stereo-ssc.nascom.nasa.gov/data/beacon//"
                      "{spacecraft}/secchi/img/euvi/%Y%m%d/%Y%m%d_%H%M%S_n7eu{spacecraft_letter}.fts")

    # instantiate each scraper
    stereo = Scraper(stereo_pattern, spacecraft=spacecraft, spacecraft_letter=spacecraft_letter)

    stereo_files = stereo.filelist(timerange)
    stereo_files.sort()

    if earliest:
        if len(stereo_files)>0:
            return stereo_files[0]
        else:
            return None
    else:
        return stereo_files

t1 = parse_time("2012-01-01").datetime
t2 = parse_time("2012-12-31").datetime

date_list = [t1]

while t2>=t1:
    t1 = t1 + relativedelta(days=1)
    date_list.append(t1)


def get_per_day(date):
	t1 = parse_time(date).datetime
	t2 = t1 + relativedelta(hours=23)
	fa = get_euvi_beacon_data(TimeRange(t1, t2), spacecraft="ahead")
	fb = get_euvi_beacon_data(TimeRange(t1, t2), spacecraft="behind")


	return [len(fa), len(fb)]
