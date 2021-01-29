import matplotlib.pyplot as plt 
import sunpy.map 
from sunpy.time import parse_time 
import numpy as np 
from dateutil.relativedelta import relativedelta
import glob

t1 = parse_time("2012-01-01").datetime
t2 = parse_time("2012-12-31").datetime

date_list = [t1]

while t2>=t1:
    t1 = t1 + relativedelta(days=1)
    date_list.append(t1)

def get_pos(date):
	date = parse_time(date)
	
	aia_file = glob.glob(date.strftime("./aia_files/*%Y_%m_%d*lev1*"))
	stereo_a_file = glob.glob(date.strftime("./stereo_files/*%Y%m%d*eua*"))
	stereo_b_file = glob.glob(date.strftime("./stereo_files/*%Y%m%d*eub*")) 

	if len(aia_file)>0 and len(stereo_a_file)>0 and len(stereo_b_file)>0:
		try:
			aia_map = sunpy.map.Map(aia_file)
			stereo_a_map = sunpy.map.Map(stereo_a_file)
			stereo_b_map = sunpy.map.Map(stereo_b_file)
		except:
			return 

	aia_coord = aia_map.observer_coordinate
	stereo_a_coord = stereo_a_map.observer_coordinate
	stereo_b_coord = stereo_b_map.observer_coordinate

	return [aia_coord, stereo_a_coord, stereo_b_coord]