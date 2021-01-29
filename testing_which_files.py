from sunpy.net import Fido, attrs as a 
import matplotlib.pyplot as plt 
import numpy as np 
import glob
from dateutil.relativedelta import relativedelta
from sunpy.time import parse_time
import re 
import datetime

files_a = glob.glob("./final_ahead/*.fts")
files_a.sort()
files_b = glob.glob("./final_behind/*.fts")
files_b.sort()


def find_file_times(file_name):
    if isinstance(file_name, (list, np.array)):
        t_list = []
        for i in range(len(file_name)):
            t = re.search('\d{8}_\d{6}', file_name[i]).group()
            tt = datetime.datetime.strptime(t, '%Y%m%d_%H%M%S')
            t_list.append(tt)
        return t_list
    else:

        t = re.search('\d{8}_\d{6}', file_name).group()
        tt = datetime.datetime.strptime(t, '%Y%m%d_%H%M%S')
        return tt


times_a = find_file_times(files_a)
times_b = find_file_times(files_b)