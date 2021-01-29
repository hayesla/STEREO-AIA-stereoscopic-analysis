import sunpy.map 
from sunpy.time import parse_time
from sunpy.map.sources.source_type import source_stretch
from sunpy.coordinates import frames, get_body_heliographic_stonyhurst
import matplotlib.pyplot as plt 
import numpy as np 
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import glob
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from dateutil.relativedelta import relativedelta
import pylab
from matplotlib import colors
import os
import subprocess
import seaborn as sns 
sns.set_style("dark")

t1 = parse_time("2011-01-01").datetime
t2 = parse_time("2014-09-29").datetime

date_list = [t1]

while t2>=t1:
    t1 = t1 + relativedelta(days=1)
    date_list.append(t1)


def get_car_map(date, prep=True):
    
    aia_file = glob.glob(date.strftime("./aia_files/*%Y_%m_%d*lev1*"))
    
    if prep:
        stereo_a_file = glob.glob(date.strftime("./prepped_a/*%Y%m%d*euA*"))
        stereo_b_file = glob.glob(date.strftime("./prepped_b/*%Y%m%d*euB*"))
    else:
        stereo_a_file = glob.glob(date.strftime("./final_ahead/*%Y%m%d*euA*"))
        stereo_b_file = glob.glob(date.strftime("./final_behind/*%Y%m%d*euB*")) 


    maps = sunpy.map.Map([aia_file, stereo_a_file, stereo_b_file])
    
    if maps[0].detector == "AIA":
        maps[0] = sunpy.map.Map(maps[0].data/maps[0].exposure_time.value/1.9, maps[0].meta)

    maps = [m.resample((512, 512)*u.pix) for m in maps]
    
    shape_out = (180, 360)  

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


    filename = "./final_maps_fixed/mosaic_{:s}.fits".format(date.strftime("%Y%m%d"))
    if not os.path.exists(filename):
        outmap.save(filename)
    

def do_all():
    errors = []
    for i in range(len(date_list)):
        print(i)
        try:
            get_car_map(date_list[i])
        except:
            print("error")
            errors.append(i)

## there were no errors from this - huzzah!

filename = "./final_maps_fixed/mosaic_{:s}.fits"
def plot_tests(i):
    
    file = filename.format(date_list[i].strftime("%Y%m%d"))
    outmap = sunpy.map.Map(file)

    plot_settings = {"cmap": euvi_map, 
                     "norm": ImageNormalize(stretch=source_stretch(outmap.meta, PowerStretch(0.20)), clip=False)}
    outmap.plot_settings = plot_settings
    cmap = outmap.cmap
    cmap.set_bad(color="k")
    outmap.nickname = 'AIA + EUVI/A + EUVI/B'
    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(projection=outmap)
    outmap.plot(vmin=-10, vmax=2000)
    plt.savefig("./tests_final_3_plots/mos_map_{:04d}.png".format(i), dpi=200)

    plt.close()


def make_movie():
    subprocess.call(['ffmpeg','-r', '20' ,'-f', 'image2', '-s', 
                    '1920x1080', '-i', "./tests_final/mos_map_%04d.png", 
                    '-vcodec', 'libx264', '-crf', '25',  '-pix_fmt', 'yuv420p', 'tests_mos_final.mp4'])


fake_data = np.zeros((512, 512))
my_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime="2017-08-01",
                    observer='earth', frame=frames.Helioprojective)
fake_header = sunpy.map.make_fitswcs_header(fake_data, my_coord, scale=[6,6]*u.arcsec/u.pixel)
fake_map = sunpy.map.Map(fake_data, fake_header)

def make_submap(mapy):
    bl = SkyCoord(-1200*u.arcsec, -1200*u.arcsec, frame=mapy.coordinate_frame)
    tr = SkyCoord(1200*u.arcsec, 1200*u.arcsec, frame=mapy.coordinate_frame)
    return mapy.submap(bl, top_right=tr)

def plot_data_final(i, date="2011-11-07"):

    date = parse_time(date_list[i])



    # get data 
    file = filename.format(date.strftime("%Y%m%d"))
    outmap = sunpy.map.Map(file)


    # get the files
    aia_file = glob.glob(date.strftime("./aia_files/*%Y_%m_%d*lev1*"))

    stereo_a_file = glob.glob(date.strftime("./prepped_a/*%Y%m%d*euA*"))
    stereo_b_file = glob.glob(date.strftime("./prepped_b/*%Y%m%d*euB*"))

    def makey_mapy(ff):
        if len(ff)>0:
            mapy = sunpy.map.Map(ff)
            return mapy.resample((512, 512)*u.pix)
        else:
            return fake_map

    map_aia = makey_mapy(aia_file)
    map_stereo_a = makey_mapy(stereo_a_file)
    map_stereo_b = makey_mapy(stereo_b_file)



    if map_aia.exposure_time>0:
        map_aia = sunpy.map.Map(map_aia.data/map_aia.exposure_time.value, map_aia.meta)



    # plotting 
    plot_settings = {"cmap": euvi_map, 
                     "norm": ImageNormalize(stretch=source_stretch(outmap.meta, PowerStretch(0.2)), clip=False)}

    outmap.plot_settings = plot_settings
    outmap.nickname = 'AIA + EUVI/A + EUVI/B'
    cmap = outmap.cmap
    cmap.set_bad(color="k")

    """
    Set up plot


    """

    map_aa = map_stereo_a.rotate()
    map_bb = map_stereo_b.rotate()
    map_aia = map_aia.rotate()


    def get_plot_lims(mapy):


        xlims_world = [-1200, 1200]*u.arcsec
        ylims_world = [-1200, 1200]*u.arcsec


        world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=mapy.coordinate_frame)
        pixel_coords = mapy.world_to_pixel(world_coords)

        # we can then pull out the x and y values of these limits.
        xlims_pixel = pixel_coords.x.value
        ylims_pixel = pixel_coords.y.value
        return xlims_pixel, ylims_pixel

    xlims_aa, ylims_aa = get_plot_lims(map_aa)
    xlims_bb, ylims_bb = get_plot_lims(map_bb)
    xlims_aia, ylims_aia = get_plot_lims(map_aia)

    fig = plt.figure(figsize=(12, 8))

    #left, bottom, width, height
    w_map = 0.28
    map_bottom = 0.68

    # map stereo_a
    ax1 = pylab.axes([0.05, map_bottom, w_map, w_map], projection=map_bb)
    # map aia
    ax2 = pylab.axes([0.25, map_bottom, w_map, w_map], projection=map_aia)
    # map stereo_b
    ax3 = pylab.axes([0.45 , map_bottom, w_map, w_map], projection=map_aa)
    # map car map
    ax4 = pylab.axes([0.05, 0.1, 0.7, 0.5], projection=outmap)
    # polar plot
    ax5 = pylab.axes([0.7, 0.6, 0.32, 0.32], projection='polar')


    # plot STEREO-B
    
    map_bb.plot(axes=ax1, cmap=euvi_map, vmin=-10, vmax=2000)
    ax1.coords.grid(lw=0.3)
    map_bb.draw_limb(lw=0.5)
    # ax1.plot_coord(SkyCoord(0*u.arcsec, 0*u.arcsec, frame=map_bb.coordinate_frame), color="r", marker="x")
    ax1.set_ylabel("HPC (Y)")
    ax1.set_xlabel("HPC (X)")
    ax1.set_title("STEREO-B/EUVI")
    ax1.set_xlim(xlims_bb)
    ax1.set_ylim(ylims_bb)
    # plot AIA 
    
    map_aia.plot(axes=ax2, vmin=0, vmax=8000)
    ax2.coords.grid(lw=0.3)
    map_aia.draw_limb(lw=0.5)
    # ax2.plot_coord(SkyCoord(0*u.arcsec, 0*u.arcsec, frame=map_aia.coordinate_frame), color="orange", marker="x")
    ax2.set_ylabel(" ")
    ax2.set_xlabel(" ")
    ax2.tick_params(axis="y", labelleft=False)
    ax2.set_title("SDO/AIA")
    ax2.set_xlim(xlims_aia)
    ax2.set_ylim(ylims_aia)
    # plot STEREO-A
    
    map_aa.plot(axes=ax3, cmap=euvi_map, vmin=-10, vmax=2000)
    ax3.coords.grid(lw=0.3)
    map_aa.draw_limb(lw=0.5)
    # ax3.plot_coord(SkyCoord(0*u.arcsec, 0*u.arcsec, frame=map_aa.coordinate_frame), color="b", marker="x")
    # ax3.coords[0].set_ticks(spacing=100*u.arcsec, color='k')
    # ax3.coords[1].set_ticks(spacing=100*u.arcsec, color='k')
    ax3.set_ylabel(" ")
    ax3.set_xlabel(" ")
    ax3.tick_params(axis="y", labelleft=False)
    ax3.set_title("STEREO-A/EUVI")
    ax3.set_xlim(xlims_aa)
    ax3.set_ylim(ylims_aa)

    # plot the HGS map
    outmap.plot(axes=ax4, vmin=-10, vmax=2000)
    ax4.set_title("AIA 193 $\mathrm{\AA}$ + EUVI/A + EUVI/B 195 $\mathrm{\AA}$ "+ str(outmap.date.datetime.strftime("%Y-%m-%d"))+" ")
    ax4.set_ylabel("Heliographic Latitude (deg)")
    ax4.set_xlabel("Heliographic Longitude (deg)")
    ax4.coords[0].set_ticks(spacing=60*u.deg, color='k')
    ax4.coords[1].set_ticks(spacing=30*u.deg, color='k')
    # ax4.plot_coord(SkyCoord(0*u.arcsec, 0*u.arcsec, frame=map_bb.coordinate_frame), color="r", marker="x")
    # ax4.plot_coord(SkyCoord(0*u.arcsec, 0*u.arcsec, frame=map_aia.coordinate_frame), color="orange", marker="x")
    # ax4.plot_coord(SkyCoord(0*u.arcsec, 0*u.arcsec, frame=map_aa.coordinate_frame), color="b", marker="x")


    for axes in (ax1, ax2, ax3, ax4):
        axes.tick_params(which="both", direction="in", color="w")

    # plot polar plot of positions
    r_unit = u.AU
    circle = plt.Circle((0.0, 0.0), (20*u.Rsun).to_value(r_unit),
                        transform=ax5.transProjectionAffine + ax5.transAxes, color="gold",
                        alpha=1, label="Sun")


    sun = get_body_heliographic_stonyhurst("sun", date)
    ax5.add_artist(circle)
    ax5.text(0, 1.15, "Earth", color="tab:blue")

    ax5.text(180*u.deg.to(u.rad), 0.1, "Sun", color="gold")

    ax5.plot(0, 1., "o", color="tab:blue")
    ax5.plot(map_bb.observer_coordinate.lon.to('rad'), 
             map_bb.observer_coordinate.radius.to(r_unit), 'x', label=map_bb.observatory, color="r")

    ax5.text(map_bb.observer_coordinate.lon.to_value('rad') + 10*u.deg.to(u.rad), 
             map_bb.observer_coordinate.radius.to_value(r_unit), "B", color="r")

    ax5.plot(map_aa.observer_coordinate.lon.to('rad'), 
             map_aa.observer_coordinate.radius.to(r_unit), 'x', label=map_aa.observatory, color="b")

    ax5.text(map_aa.observer_coordinate.lon.to_value('rad') - 10*u.deg.to(u.rad), 
             map_aa.observer_coordinate.radius.to_value(r_unit), "A", color="b")


    ax5.plot(map_aia.observer_coordinate.lon.to('rad'), 
             map_aia.observer_coordinate.radius.to(r_unit), 'x', label=map_aia.observatory, color="orange")

    ax5.text(map_aia.observer_coordinate.lon.to_value('rad') + 2.5*u.deg.to(u.rad), 
             map_aia.observer_coordinate.radius.to_value(r_unit), "SDO", color="orange")

    ax5.text(map_aia.observer_coordinate.lon.to_value('rad') + 2.5*u.deg.to(u.rad), 
             map_aia.observer_coordinate.radius.to_value(r_unit), "SDO", color="orange")

    ax5.plot([map_bb.observer_coordinate.lon.to_value('rad'), sun.lon.to_value('rad')], 
            [map_bb.observer_coordinate.radius.to_value(r_unit), sun.radius.to_value(r_unit)], 
            ls="dashed", color="k", lw=0.3)

    ax5.plot([map_aa.observer_coordinate.lon.to_value('rad'), sun.lon.to_value('rad')], 
            [map_aa.observer_coordinate.radius.to_value(r_unit), sun.radius.to_value(r_unit)],
            ls="dashed", color="k", lw=0.3)

    ax5.plot([map_aia.observer_coordinate.lon.to_value('rad'), sun.lon.to_value('rad')], 
            [map_aia.observer_coordinate.radius.to_value(r_unit), sun.radius.to_value(r_unit)],
            ls="dashed", color="k", lw=0.3)

    ax5.set_theta_zero_location("S")
    ax5.set_rlim(0, 1.3)
    ax5.set_rlabel_position(45.5)
    ax5.legend(bbox_to_anchor=(1.1, -0.3), loc="lower right")

    plt.savefig("./tests_final/mos_map_{:04d}.png".format(i), dpi=200)
    plt.close()

def test_aia(date, **kwargs):
    plt.close()
    aia_file = glob.glob(date.strftime("./aia_files/*%Y_%m_%d*lev1*"))
    map_aia = sunpy.map.Map(aia_file)

    map_aia = sunpy.map.Map(map_aia.data/map_aia.exposure_time.value, map_aia.meta)
    map_aia.plot(**kwargs)
    plt.colorbar()
    print(map_aia.data.min(), map_aia.data.max())



def plot_polar(i):


    date = parse_time(date_list[i])

    aia_file = glob.glob(date.strftime("./aia_files/*%Y_%m_%d*lev1*"))
    stereo_a_file = glob.glob(date.strftime("./prepped_a/*%Y%m%d*euA*"))
    stereo_b_file = glob.glob(date.strftime("./prepped_b/*%Y%m%d*euB*"))
    
    def makey_mapy(ff):
        if len(ff)>0:
            mapy = sunpy.map.Map(ff)
            return mapy.resample((512, 512)*u.pix)
        else:
            return fake_map

    map_aia = makey_mapy(aia_file)
    map_aa = makey_mapy(stereo_a_file)
    map_bb = makey_mapy(stereo_b_file)


    sun = get_body_heliographic_stonyhurst("sun", date)
    plt.figure(figsize=(8, 8))
    ax = plt.subplot(projection='polar')
    r_unit = u.AU

    ax.text(0, 1.05, "Earth", color="tab:blue")

    ax.text(0.2, 0, "Sun", color="yellow")

    ax.plot(0, 1., "o", color="tab:blue")


    ax.plot(map_bb.observer_coordinate.lon.to('rad'), 
             map_bb.observer_coordinate.radius.to(r_unit), 'x', label=map_bb.observatory, color="r")
    ax.plot(map_aa.observer_coordinate.lon.to('rad'), 
             map_aa.observer_coordinate.radius.to(r_unit), 'x', label=map_aa.observatory, color="b")
    ax.plot(map_aia.observer_coordinate.lon.to('rad'), 
             map_aia.observer_coordinate.radius.to(r_unit), 'x', label=map_aia.observatory, color="orange")

    ax.text(map_aia.observer_coordinate.lon.to_value('rad')+ 2.5*u.deg.to(u.rad), 
             map_aia.observer_coordinate.radius.to_value(r_unit), "SDO", color="orange")

    ax.text(map_aa.observer_coordinate.lon.to_value('rad')- 5*u.deg.to(u.rad), 
             map_aa.observer_coordinate.radius.to_value(r_unit), "A", color="b")

    ax.text(map_bb.observer_coordinate.lon.to_value('rad')+ 5*u.deg.to(u.rad), 
             map_bb.observer_coordinate.radius.to_value(r_unit), "B", color="r")

    ax.plot([map_bb.observer_coordinate.lon.to_value('rad'), sun.lon.to_value('rad')], 
            [map_bb.observer_coordinate.radius.to_value(r_unit), sun.radius.to_value(r_unit)], 
            ls="dashed", color="k", lw=0.3)

    ax.plot([map_aa.observer_coordinate.lon.to_value('rad'), sun.lon.to_value('rad')], 
            [map_aa.observer_coordinate.radius.to_value(r_unit), sun.radius.to_value(r_unit)],
            ls="dashed", color="k", lw=0.3)

    ax.plot([map_aia.observer_coordinate.lon.to_value('rad'), sun.lon.to_value('rad')], 
            [map_aia.observer_coordinate.radius.to_value(r_unit), sun.radius.to_value(r_unit)],
            ls="dashed", color="k", lw=0.3)
    ax.set_theta_zero_location("S")


    circle = plt.Circle((0.0, 0.0), (10*u.Rsun).to_value(r_unit),
                        transform=ax.transProjectionAffine + ax.transAxes, color="yellow",
                        alpha=1, label="Sun")

    ax.add_artist(circle)
    

    plt.show()

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