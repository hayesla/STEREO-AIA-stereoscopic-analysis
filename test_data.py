import sunpy.map 
import matplotlib.pyplot as plt 
import numpy as np
import glob
from astropy import units as u 
import subprocess 

aia_files = glob.glob("./aia_files/*.fits")
aia_files.sort()

def test_aia():
	for i in range(534, len(aia_files)):
		print(i)
		aia_map = sunpy.map.Map(aia_files[i])
		fig = plt.figure()
		ax  = fig.add_subplot(projection=aia_map)
		aia_map.plot(clip_interval=(1., 99.99)*u.percent)

		plt.savefig("./test_plots/aia_map_{:04d}.png".format(i))
		plt.close()

	subprocess.call(['ffmpeg','-r', '10' ,'-f', 'image2', '-s', 
					'1920x1080', '-i', "./test_plots/aia_map_%04d.png", 
					'-vcodec', 'libx264', '-crf', '25',  '-pix_fmt', 'yuv420p', 'aia_test_mov.mp4'])



fake_map = sunpy.map.Map(aia_files[0])

stereo_a_files = glob.glob("./stereo_files/*eua.fts")
stereo_a_files.sort()



stereo_b_files = glob.glob("./stereo_files/*eub.fts")
stereo_b_files.sort()


def test_stereo_a():
	errors = []
	for i in range(0, len(stereo_a_files)):
		try:
			print(i)
			stereo_map = sunpy.map.Map(stereo_a_files[i])
			fig = plt.figure()
			ax  = fig.add_subplot(projection=stereo_map)
			stereo_map.plot(clip_interval=(1., 99.99)*u.percent)

			plt.savefig("./test_plots/stereo_a_{:04d}.png".format(i))
			plt.close()
		except:
			errors.append(i)

			fig = plt.figure()
			ax  = fig.add_subplot(projection=fake_map)
			fake_map.plot(clip_interval=(1., 99.99)*u.percent)

			plt.savefig("./test_plots/stereo_a_{:04d}.png".format(i))
			plt.close()

	subprocess.call(['ffmpeg','-r', '10' ,'-f', 'image2', '-s', 
					'1920x1080', '-i', "./test_plots/stereo_a_%04d.png", 
					'-vcodec', 'libx264', '-crf', '25',  '-pix_fmt', 'yuv420p', 'stereo_a_test_mov.mp4'])


def test_stereo_b():
	errors = []
	for i in range(0, len(stereo_b_files)):
		try:
			print(i)
			stereo_map = sunpy.map.Map(stereo_b_files[i])
			fig = plt.figure()
			ax  = fig.add_subplot(projection=stereo_map)
			stereo_map.plot(clip_interval=(1., 99.99)*u.percent)

			plt.savefig("./test_plots/stereo_b_{:04d}.png".format(i))
			plt.close()
		except:
			errors.append(i)

			fig = plt.figure()
			ax  = fig.add_subplot(projection=fake_map)
			fake_map.plot(clip_interval=(1., 99.99)*u.percent)

			plt.savefig("./test_plots/stereo_b_{:04d}.png".format(i))
			plt.close()

	subprocess.call(['ffmpeg','-r', '10' ,'-f', 'image2', '-s', 
					'1920x1080', '-i', "./test_plots/stereo_b_%04d.png", 
					'-vcodec', 'libx264', '-crf', '25',  '-pix_fmt', 'yuv420p', 'stereo_b_test_mov.mp4'])
	return errors