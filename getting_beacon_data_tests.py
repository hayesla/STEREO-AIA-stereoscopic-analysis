import urllib
from urllib.error import HTTPError, URLError
from bs4 import BeautifulSoup , SoupStrainer
import re
import datetime
from dateutil.relativedelta import relativedelta
from sunpy.time import parse_time
import numpy as np

savedir_a = "/Users/laurahayes/spaceweather_stuff/fun_stuff/final_ahead/"
savedir_b = "/Users/laurahayes/spaceweather_stuff/fun_stuff/final_behind/"

t1 = parse_time("2011-01-01").datetime
t2 = parse_time("2014-10-01").datetime
date_list = [t1]

from parfive import Downloader

while t2>=t1:
    t1 = t1 + relativedelta(days=1)
    date_list.append(t1)

base_url = "https://www.solarmonitor.org/data/"

url_tests_b = [d.strftime(base_url+"%Y/%m/%d/fits/strb/") for d in date_list]
url_tests_a = [d.strftime(base_url+"%Y/%m/%d/fits/stra/") for d in date_list]

beacon_url_a = "https://stereo-ssc.nascom.nasa.gov/data/beacon//ahead/secchi/img/euvi/"
beacon_url_b = "https://stereo-ssc.nascom.nasa.gov/data/beacon//behind/secchi/img/euvi/"


beacon_test_a = [d.strftime(beacon_url_a+"%Y%m%d/") for d in date_list]
beacon_test_b = [d.strftime(beacon_url_b+"%Y%m%d/") for d in date_list]


def list_path_files(url_path, file_name):
    test = urllib.request.urlopen(url_path)
    soup = BeautifulSoup(test, features="lxml")

    fits_links = []
    for link in soup.findAll('a'):
        if link.get('href') is not None and link.get('href').find(file_name) != -1:
            fits_links.append(url_path + link.get('href').split('/')[-1]) 

    return fits_links

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

def get_urls_hours(beacon_tests):
    files_beacon = []
    for u in beacon_tests:
        print(u)
        try:
            urls = list_path_files(u, "fts")
            files_beacon.append(urls)
        except HTTPError:
            print("error")
            files_beacon.append([])

    datess = [find_file_times(x) for x in files_beacon]
    hours = []
    for d in datess:
        try:
            hours.append(parse_time(d).strftime("%H").astype(int))
        except:
            hours.append([])

    date_list_12 = [d + relativedelta(hours=12) for d in date_list]

    testy = []
    final_urls = []
    for i in range(len(datess)):
        if len(datess[i])>0:
            ind = np.argmin(np.abs(np.array(datess[i]) - date_list_12[i]))
            testy.append(datess[i][ind])
            final_urls.append(files_beacon[i][ind])
        else:
            testy.append([])
            final_urls.append([])
    return files_beacon, final_urls, testy

# all_ahead, final_ahead, time_ahead = get_urls_hours(beacon_test_a)
# all_behind, final_behind, time_behind = get_urls_hours(beacon_test_b)

def save_files():
    lala = {"file_list":files_beacon_ahead, "final_urls":final_urls, "dates_final":testy}
    lala2 = {"file_list":all_behind, "final_urls":final_behind, "dates_final":time_behind}
    with open("ahead_files.pkl", "wb") as handle:
        pickle.dump(lala, handle)

    with open("behind_files.pkl", "wb") as handle:
        pickle.dump(lala2, handle)


def download_files(file_list, savedir):
    dl = Downloader()
    for f in file_list:
        dl.enqueue_file(f, path=savedir)
    files = dl.download()

def download_files2(file_list, savedir):
    for f in file_list:
        if isinstance(f, str):
            fname = savedir + f.split("/")[-1]
            urllib.request.urlretrieve(f, fname)
        else:
            print(f)


# filess = []
# for u in url_tests_a:
#   print(u)
#   try:
#       urls = list_path_files(u, "fts.gz")
#       filess.append(urls)
#   except HTTPError:
#       print("error")
#       filess.append([])