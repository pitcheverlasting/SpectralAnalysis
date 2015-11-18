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
	[ax.semilogx(fs/freq, sqrt(coh[i, :]), color=col[i], linewidth=2.5, label=drivername[i]) for i in xrange(0, nvar)]
	# ax.xaxis.set_visible(False)
	# plt.xlim([styr, edyr]) # so that no space before and after the time series
	# ax.set_xticks(years[2::20])  # specify where you put the ticks
	# ax.set_xticklabels(years[2::20], fontsize=15)  # adjust the ticks font size
	# ax.set_ylabel("Z-score", fontsize=24)
	# ax.tick_params(axis='y', labelsize=18)
	ymin, ymax = ax.get_ylim()
	ax.text(1050, ymax-0.2*(ymax-ymin), "%s" % title, fontsize=22)
	ax.plot([365, 365], [0, 1], 'r--', linewidth=3.5)
	ax.plot([120, 120], [0, 1], 'k--', linewidth=3.5)
	ax.plot([30, 30], [0, 1], 'b--', linewidth=3.5)
	plt.grid(True)



	return