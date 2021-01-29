import matplotlib.pyplot as plt
import numpy as np
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
import sunpy.sun
from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
import datetime
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from sunpy.map.sources.source_type import source_stretch
from dateutil.relativedelta import relativedelta
import glob 
import matplotlib.colors as colors

t1 = parse_time("2012-01-01").datetime
t2 = parse_time("2012-12-31").datetime

date_list = [t1]

while t2>=t1:
    t1 = t1 + relativedelta(days=1)
    date_list.append(t1)


def plot_data(i, date="2011-11-07"):

    print(date)
    date = parse_time(date)

    # get the files
    aia_file = glob.glob(date.strftime("./aia_files/*%Y_%m_%d*lev1*"))
    stereo_a_file = glob.glob(date.strftime("./stereo_files/*%Y%m%d*eua*"))
    stereo_b_file = glob.glob(date.strftime("./stereo_files/*%Y%m%d*eub*")) 

    files = [aia_file, stereo_a_file, stereo_b_file]
    
    # make maps
    maps = sunpy.map.Map(sorted(files))
    maps = [m.resample((1024, 1024)*u.pix) for m in maps]

    shape_out = (180, 360)  # This is set deliberately low to reduce memory consumption

    header = sunpy.map.make_fitswcs_header(shape_out,
                                           SkyCoord(0, 0, unit=u.deg,
                                                    frame="heliographic_stonyhurst",
                                                    obstime=maps[0].date),
                                           scale=[180 / shape_out[0],
                                                  360 / shape_out[1]] * u.deg / u.pix,
                                           wavelength=int(maps[0].meta['wavelnth']) * u.AA,
                                           projection_code="CAR")
    out_wcs = WCS(header)

    coordinates = tuple(map(sunpy.map.all_coordinates_from_map, maps))
    weights = [coord.transform_to("heliocentric").z.value for coord in coordinates]

    weights = [(w / np.nanmax(w)) ** 3 for w in weights]
    for w in weights:
        w[np.isnan(w)] = 0

    array, _ = reproject_and_coadd(maps, out_wcs, shape_out,
                                   input_weights=weights,
                                   reproject_function=reproject_interp,
                                   match_background=True,
                                   background_reference=0)


    outmap = sunpy.map.Map((array, header))

    # plotting 
    plot_settings = {"cmap": "sohoeit195", 
                 "norm": ImageNormalize(stretch=source_stretch(outmap.meta, PowerStretch(0.25)), clip=False)}

    outmap.plot_settings = plot_settings
    outmap.nickname = 'AIA + EUVI/A + EUVI/B'

    cmap = outmap.cmap
    cmap.set_bad(color="k")

    fig = plt.figure()
    plt.figure(figsize=(10, 5))
    ax = plt.subplot(projection=out_wcs)
    im = outmap.plot(vmin=400)

    outmap.save("./final_maps/mosaic_{:s}.fits".format(date.strftime("%Y%m%d")))


    plt.savefig(date.strftime("./plots/aia_stereo_{:03d}.png".format(i)))
    plt.close()

fake_plot = sunpy.map.Map("mosaic_20120101.fits")
plot_settings = {"cmap": "sohoeit195", 
             "norm": ImageNormalize(stretch=source_stretch(fake_plot.meta, PowerStretch(0.25)), clip=False)}
fake_plot.plot_settings = plot_settings

def do_all():
    errors = []
    
    for i in range(132, len(date_list)):
        try:
            plot_data(i, date=date_list[i])
        except:
            fig = plt.figure()
            plt.figure(figsize=(10, 5))
            ax = plt.subplot(projection=fake_plot)
            im = fake_plot.plot(vmin=400)
            plt.savefig(date.strftime("./plots/aia_stereo_{:03d}.png".format(i)))
            plt.close()
            errors.append(d)
    return errors

def make_movie():
	subprocess.call(['ffmpeg','-r', '10' ,'-f', 'image2', '-s', 
					'1920x1080', '-i', "./plots/aia_stereo_%03d.png", 
					'-vcodec', 'libx264', '-crf', '25',  '-pix_fmt', 'yuv420p', 'stereo_aia_2012.mp4'])

def do_from_saved():
    files = glob.glob("./final_maps/*.fits")
    files.sort()
    for i in range(len(files)):
        print(i)
        outmap = sunpy.map.Map(files[i])
        plot_settings = {"cmap": euvi_map, 
                         "norm": ImageNormalize(stretch=source_stretch(outmap.meta, PowerStretch(0.25)), clip=False),
                         "interpolation": "bilinear"}

        outmap.plot_settings = plot_settings
        cmap = outmap.cmap
        cmap.set_bad(color="k")

        fig = plt.figure(figsize=(10, 5))
        ax = plt.subplot(projection=outmap)
        im = outmap.plot(clip_interval=(1, 99.99)*u.percent)
        plt.savefig("./plots/aia_stereo_{:03d}.png".format(i))
        plt.close()
    
    subprocess.call(['ffmpeg','-r', '10' ,'-f', 'image2', '-s', 
                    '1920x1080', '-i', "./plots/aia_stereo_%03d.png", 
                    '-vcodec', 'libx264', '-crf', '25',  '-pix_fmt', 'yuv420p', 'stereo_aia_2012_new_cmap.mp4'])


def cmap_from_rgb_file(name, fname):
    """
    Create a colormap from a RGB .csv file.
    The .csv file must have 3  equal-length columns of integer data, with values
    between 0 and 255, which are the red, green, and blue values for the colormap.
    Parameters
    ----------
    name : str
        Name of the colormap.
    fname : str
        Filename of data file. Relative to the sunpy colormap data directory.
    Returns
    -------
    cmap : matplotlib.colors.LinearSegmentedColormap
    """
    data = np.loadtxt(fname, delimiter=',')
    if data.shape[1] != 3:
        raise RuntimeError(f'RGB data files must have 3 columns (got {data.shape[1]})')
    return _cmap_from_rgb(data[:, 0], data[:, 1], data[:, 2], name)


def _cmap_from_rgb(r, g, b, name):
    cdict = create_cdict(r, g, b)
    return colors.LinearSegmentedColormap(name, cdict)


def create_cdict(r, g, b):
    """
    Create the color tuples in the correct format.
    """
    i = np.linspace(0, 1, r.size)

    cdict = {name: list(zip(i, el / 255.0, el / 255.0))
             for el, name in [(r, 'red'), (g, 'green'), (b, 'blue')]}
    return cdict


euvi_map = cmap_from_rgb_file("euvi195", "/Users/laurahayes/sunpy_dev/euvi_195.csv")



