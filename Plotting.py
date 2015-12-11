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
	drivername = ['Tair','Rs','Wind','Humidity','Cloudiness', 'Tsurf', 'VPD', 'Rain']
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

def CoherenceBasinPlotTogether(ax, coh, fs, freq, basinnames):

	nvar = coh.shape[0]
	# col = ['LightCoral', 'Orange', 'MediumSpringGreen', 'Orchid', 'SkyBlue', 'Crimson', 'Purple', 'Aquamarine', 'red', 'black']
	cm_sub = linspace(0, 1, 10)
	col = [plt.cm.gist_rainbow(x) for x in cm_sub]
	# linestyle = ['-o','-s','-*','-^','-','-o','-s','-*']
	[ax.semilogx(fs/freq, coh[i, :], color=col[i], linewidth=5, alpha=0.2, label=basinnames[i]) for i in xrange(0, nvar)]
	# ax.xaxis.set_visible(False)
	plt.xlim([1, 2000]) # so that no space before and after the time series
	# ax.set_xticks(years[2::20])  # specify where you put the ticks
	# ax.set_xticklabels(years[2::20], fontsize=15)  # adjust the ticks font size
	# ax.set_ylabel("Z-score", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	# ax.set_title("%s" % title, fontsize=18, y=1.02)
	# ymin, ymax = ax.get_ylim()
	# ax.text(1.5, ymax-0.15*(ymax-ymin), "%s" % title, fontsize=20)
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	ax.plot([90, 90], [0, 1], 'k--', linewidth=3.5)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.grid(True)

	return

def Coherence_obs_mod_Plot(ax, coh, fs, freq, title):

	# linestyle = ['-o','-s','-*','-^','-','-o','-s','-*']
	ax.semilogx(fs/freq, coh, color='grey', linewidth=3)
	# ax.xaxis.set_visible(False)
	plt.xlim([1, 2000]) # so that no space before and after the time series
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
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

def Mapshow(data, lons, lats, size, alpha, min, max, cmp, title, legend, unit, figdir, filename):

	m = Basemap(projection='merc', llcrnrlon=71, llcrnrlat=3, urcrnrlon=139, urcrnrlat=55, lat_ts=20)
	# draw boundaries using shapefile with countries boundaries
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=2)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=1)
	m.drawparallels(np.arange(10, 60, 20), labels=[1, 0, 0, 0], fontsize=19)  # only left ytick
	m.drawmeridians(np.arange(90, 140, 30), labels=[0, 0, 0, 1], fontsize=19)  # only bottom xtick
	# draw data
	if type(data) is str:
		im = m.scatter(lons, lats, size, marker='o', c=data, latlon=True, label=legend)
	else:
		im = m.scatter(lons, lats, size, marker='o', c=data, edgecolors='none', alpha=alpha, vmin=min, vmax=max, latlon=True, cmap=cmp, label=legend)
		cb = m.colorbar(im, pad='3%')
		cb.ax.tick_params(labelsize=20)
	# plotdata = m.transform_scalar(data, lons, lats, nx, ny)
	# im = m.imshow(plotdata, vmin=min, vmax=max, cmap=cmp)
	plt.title(title, fontsize=16, y=1.02)
	# plt.legend(loc=3, prop={'size': 18})
	plt.xlabel(unit, fontsize=21, labelpad=19)
	# savefig('%s%s' % (figdir, filename))
	plt.show()
	# plt.clf()

def Mapshow_RGB(data, lons, lats, size, alpha, marker, title, legend, unit):

	m = Basemap(projection='merc', llcrnrlon=71, llcrnrlat=3, urcrnrlon=139, urcrnrlat=55, lat_ts=20)
	# draw boundaries using shapefile with countries boundaries
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=2)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=1)
	m.drawparallels(np.arange(10, 60, 20), labels=[1, 0, 0, 0], fontsize=19)  # only left ytick
	m.drawmeridians(np.arange(90, 140, 30), labels=[0, 0, 0, 1], fontsize=19)  # only bottom xtick
	# draw data
	im = m.scatter(lons, lats, size, marker=marker, c=data, edgecolors=data, alpha=alpha, latlon=True, label=legend)
	# cb = m.colorbar(im, pad='3%')
	# plotdata = m.transform_scalar(data, lons, lats, nx, ny)
	# im = m.imshow(plotdata, vmin=min, vmax=max, cmap=cmp)
	plt.title(title, fontsize=16, y=1.02)
	plt.xlabel(unit, labelpad=20)
	# plt.show()
	# plt.clf()

def MapShow_lambert():
	m = Basemap(llcrnrlon=82, llcrnrlat=0, urcrnrlon=140, urcrnrlat=55, projection='lcc', lat_1=20, lat_2=40, lon_0=108)
	# draw boundaries using shapefile with countries boundaries
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou1_4m/bou1_4l', 'CHN_adm1', linewidth=2)
	m.readshapefile('/home/water5/lpeng/Masks/Shapefile/china_map/bou2_4m/bou2_4l', 'CHN_adm2', linewidth=1)
	# m.etopo()
	m.drawparallels(np.arange(20, 60, 20), labels=[1, 0, 0, 0])  # only left ytick
	m.drawmeridians(np.arange(90, 140, 30), labels=[0, 0, 0, 1])  # only bottom xtick
	plt.show()

def Stacked_Bar(ax, data, width, positions, stackbase, col, var, sign, alp):

	mask = zeros(data.shape)
	if sign == "positive":
		mask[data>0.0] = 1.0
		lab = var
	elif sign == "negative":
		mask[data<0.0] = 1.0
		lab = None
	else:
		mask = ones(data.shape)
	ax.bar(positions + width / 2, data*mask, label=lab, bottom=stackbase, width=width, linewidth=0.5, color=col, alpha=alp)
	stackbase = stackbase + data * mask

	return stackbase