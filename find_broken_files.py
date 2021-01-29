import glob
import matplotlib.pyplot as plt 
import numpy as np 
from astropy.io import fits
import sunpy.map
from sunpy.time import parse_time
from astropy import units as u
from dateutil.relativedelta import relativedelta
import pickle

t1 = parse_time("2011-01-01").datetime
t2 = parse_time("2014-09-30").datetime

date_list = [t1]

while t2>=t1:
    t1 = t1 + relativedelta(days=1)
    date_list.append(t1)

filesa_dir = "./final_ahead/"
filesb_dir = "./final_behind/"

filesa_prep_dir = "./prepped_a/"
filesb_prep_dir = "./prepped_b/"
sa_files_prep = glob.glob(filesa_prep_dir + "*euA*")
sb_files_prep = glob.glob(filesb_prep_dir + "*euB*")


filesaia_dir = "./aia_files/"

aia_files = glob.glob(filesaia_dir + "*.fits")
sa_files = glob.glob(filesa_prep_dir + "*euA*")
sb_files = glob.glob(filesb_prep_dir + "*euB*")

aia_files.sort()
sa_files.sort()
sb_files.sort()

aia_files = np.array(aia_files)
sa_files = np.array(sa_files)
sb_files = np.array(sb_files)


def get_broken_stereo_files(sat="a"):
	if sat == "a":
		filesearch="./prepped_a/*%Y%m%d*euA*"
	elif sat=="b":
		filesearch="./prepped_b/*%Y%m%d*euB*"
	inds_a = []
	inds_a2 = []
	for i in range(len(date_list)-1):
	    date = date_list[i]
	    aia_file = glob.glob(date.strftime(filesearch))
	    if len(aia_file)>0:
	        aa = fits.open(aia_file[0])
	        data = np.sum(aa[0].data)
	        if data == 0:
	            inds_a.append(aia_file[0])
	    else:
	        print("error")
	        inds_a2.append(i)
	return inds_a

def get_info(sat="a"):
	if sat=="a":
		file = "ahead_files.pkl"
	elif sat=="b":
		file = "behind_files.pkl"
	ff = open(file, "rb")
	data = pickle.load(ff)
	ff.close()
	return data

info_a = get_info(sat="a")
info_b = get_info(sat="b")

files_broken_a = get_broken_stereo_files(sat="a")
files_broken_b = get_broken_stereo_files(sat="b")

data_list_a = np.array(info_a["file_list"])[inds_a]
data_list_b = np.array(info_b["file_list"])[inds_b]


broken_a_files = [b.split("/")[-1] for b in sa_files[inds_a]]
broken_b_files = [b.split("/")[-1] for b in sa_files[inds_b]]

def get_filenames_from_full(data_list_a):
	if len(data_list_a)==0:
		return []
	else:
		listy = []
		for i in range(len(data_list_a)):
			listy.append([x.split("/")[-1] for x in data_list_a[i]])
		return listy

data_file_list_a = get_filenames_from_full(data_list_a)
data_file_list_b = get_filenames_from_full(data_list_b)

def get_indices():
	inds_listy = []
	for i in range(len(broken_a_files)):
		indy = np.where(np.array(data_file_list_a[i]) == broken_a_files[i])
		inds_listy.append(indy)
	return inds_listy