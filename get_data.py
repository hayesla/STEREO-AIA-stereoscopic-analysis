from sunpy.net import Fido, attrs as a 
from astropy import units as u 
from sunpy.util.scraper import Scraper
from dateutil.relativedelta import relativedelta
from sunpy.time import parse_time

#t1 = parse_time("2011-01-01").datetime
# t1 = parse_time("2011-03-27").datetime
# t2 = parse_time("2012-01-01").datetime

# t1 = parse_time("2012-01-01").datetime

t1 = parse_time("2014-03-01").datetime
t2 = parse_time("2015-01-01").datetime


def get_aia(t1, t2):
	aia_res = Fido.search(a.Time(t1, t2), a.Instrument("AIA"), 
						  a.Wavelength(193*u.angstrom), a.vso.Sample(24*u.hour))

	Fido.fetch(aia_res, path="./aia_files")

def get_stereo(t1, t2, stereo_a=True, stereo_b=True):
	""" 
	have to do 30 days at a time

	"""
	if stereo_a:
		stereo_a_res = Fido.search(a.Time(t1, t2), a.Instrument("EUVI"), 
							  a.Wavelength(195*u.angstrom), a.Sample(24*u.hour),
							  a.vso.Source("STEREO_A"))
		Fido.fetch(stereo_a_res, path="./stereo_files")
	
	if stereo_b:
		stereo_b_res = Fido.search(a.Time(t1, t2), a.Instrument("EUVI"), 
						  a.Wavelength(195*u.angstrom), a.Sample(24*u.hour),
						  a.vso.Source("STEREO_B"))						 	
		Fido.fetch(stereo_b_res, path="./stereo_files")	

def get_all_stereo(t1, t2, stereo_a=True, stereo_b=True):
	tstart = t1

	while tstart <= t2:
		print(tstart)
		tend = tstart + relativedelta(months=1)
		print(tend)
		get_stereo(tstart, tend, stereo_a=stereo_a, stereo_b=stereo_b)
		tstart = tstart + relativedelta(months=1)






