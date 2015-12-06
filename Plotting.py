__author__ = 'lpeng'
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.basemap import Basemap, cm
import pandas as pd
import string

def CoherenceStationPlot(coh, fs, freq):

	nvar = len(coh)

	col = ['LightCoral', 'Orange', 'MediumSpringGreen', 'Orchid', 'SkyBlue', 'Crimson', 'Purple', 'Aquamarine']
	drivername = ['Tair','Rs','Wind','Humidity','Low Cloud', 'Tsurf', 'VPD', 'Rain']
	# linestyle = ['-o','-s','-*','-^','-','-o','-s','-*']
	plt.figure()
	[plt.semilogx(fs/freq, coh[i], color=col[i], linewidth=2.5) for i in xrange(0, nvar)]
	plt.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	plt.plot([120, 120], [0, 1], 'k--', linewidth=3.5)
	plt.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.show()
	plt.clf()

	return

def CoherenceBasinPlot(ax, coh, fs, freq, title):

	nvar = coh.shape[0]
	col = ['LightCoral', 'Orange', 'MediumSpringGreen', 'Orchid', 'SkyBlue', 'Crimson', 'Purple', 'Aquamarine']
	drivername = ['Tair','Rs','Wind','Humidity','Low Cloud', 'Tsurf', 'VPD', 'Rain']
	# linestyle = ['-o','-s','-*','-^','-','-o','-s','-*']
	[ax.semilogx(fs/freq, coh[i, :], color=col[i], linewidth=2.5, label=drivername[i]) for i in xrange(0, nvar)]
	# ax.xaxis.set_visible(False)
	plt.xlim([1, 2000]) # so that no space before and after the time series
	# ax.set_xticks(years[2::20])  # specify where you put the ticks
	# ax.set_xticklabels(years[2::20], fontsize=15)  # adjust the ticks font size
	# ax.set_ylabel("Z-score", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	# ax.set_title("%s" % title, fontsize=18, y=1.02)
	ymin, ymax = ax.get_ylim()
	ax.text(1.5, ymax-0.15*(ymax-ymin), "%s" % title, fontsize=20)
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	ax.plot([90, 90], [0, 1], 'k--', linewidth=3.5)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.grid(True)

	return

def PSDBasinPlot(ax, psd, fs, freq, title):

	nvar = psd.shape[0]
	col = ['LightCoral', 'Orange', 'MediumSpringGreen', 'Orchid', 'SkyBlue', 'Crimson', 'Purple', 'Aquamarine']
	drivername = ['Tair','Rs','Wind','Humidity','Low Cloud', 'Tsurf', 'VPD', 'Rain']
	# linestyle = ['-o','-s','-*','-^','-','-o','-s','-*']
	[ax.semilogy(fs/freq, psd[i, :], color=col[i], linewidth=2.5, label=drivername[i]) for i in xrange(0, nvar)]
	# ax.xaxis.set_visible(False)
	# plt.xlim([1, 2000]) # so that no space before and after the time series
	# ax.set_xticks(years[2::20])  # specify where you put the ticks
	# ax.set_xticklabels(years[2::20], fontsize=15)  # adjust the ticks font size
	# ax.set_ylabel("Z-score", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	# ax.set_title("%s" % title, fontsize=18, y=1.02)
	# ymin, ymax = ax.get_ylim()
	# ax.text(1.5, ymax-0.15*(ymax-ymin), "%s" % title, fontsize=20)
	# ax.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	# ax.plot([90, 90], [0, 1], 'k--', linewidth=3.5)
	# ax.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.grid(True)

	return

def Mapshow(data, lons, lats, min, max, cmp, tit, unit, figdir, filename):
	m = Basemap(llcrnrlon=82, llcrnrlat=0, urcrnrlon=140, urcrnrlat=55, projection='lcc', lat_1=20, lat_2=40, lon_0=108)
	# draw boundaries
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=2)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=1)
	# m.drawparallels(np.arange(25, 65, 20), labels=[1, 0, 0, 0])  # only left ytick
	# m.drawmeridians(np.arange(-120, -40, 20), labels=[0, 0, 0, 1])  # only bottom xtick

	# idx = np.where(np.isnan(data) == 0)
	# lon = lons[idx[1]]
	# lat = lats[idx[0]]
	# data = data[idx]
	im = m.scatter(lons, lats, 50, marker='o', c=data, vmin=min, vmax=max, latlon=True, cmap=cmp)
	cb = m.colorbar(im, pad='3%')
	# plotdata = m.transform_scalar(data, lons, lats, nx, ny)
	# im = m.imshow(plotdata, vmin=min, vmax=max, cmap=cmp)
	plt.title(tit)
	plt.xlabel(unit, labelpad=20)
	# savefig('%s%s' % (figdir, filename))
	plt.show()
	plt.clf()

def ShapefileShow():
	m = Basemap(llcrnrlon=82, llcrnrlat=0, urcrnrlon=140, urcrnrlat=55, projection='lcc', lat_1=20, lat_2=40, lon_0=108)
	# m = Basemap(projection='merc', llcrnrlon=70, llcrnrlat=15, urcrnrlon=140, urcrnrlat=55, lat_ts=20)
	# draw boundaries
	# m.etopo()
	# m.drawparallels(np.arange(25, 65, 20), labels=[1, 0, 0, 0])  # only left ytick
	# m.drawmeridians(np.arange(-120, -40, 20), labels=[0, 0, 0, 1])  # only bottom xtick
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=2)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=1)
	plt.show()
