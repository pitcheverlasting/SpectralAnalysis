__author__ = 'lpeng'

from pylab import *
import pickle
import scipy.io
import scipy.stats
import pandas as pd
import IO, FFT, Plotting
from PET_Library import Data
from PETREND import TrendAnalysis, TrendAttribution
# #=============================path================================
datadir = '/home/water5/lpeng/script/PROJECT/pan_evaporation_china/original'
workspace = '/home/water5/lpeng/script/PROJECT/pan_evaporation_china/workspace'
figdir = '/home/water5/Figure/pan_spectral/trend/'
##=============================variable============================
vars = ['time', 'p', 'tavg', 'tmin', 'tmax', 'ea', 'rh', 'tc', 'lc', 'wind', 'ts', 'sun', 'rain', 'pan', 'vpd',  'estavg', 'rh_test']  # time 1,2,3 year, month, day # Note that original data is multiplied by 10 from the normal unit. This process data devide 10 so that it becomes the original unit
varlongs = ['time', 'daily averaged pressure [hpa]', 'daily averaged air temperature [degree]', 'daily minimum air temperature [degree]', 'daily maximum air temperature [degree]', 'daily averaged vapor pressure [hpa]', 'daily averaged relative humidity [%]', 'daily averaged total cloud cover', 'daily averaged low cloud cover', 'daily averaged wind speed [m/s]', 'daily averaged surface temperature(at 0cm level) [degree]', 'daily averaged sunshine hour [hr]', 'daily averaged rainfall [mm]', 'daily averaged pan evaporation [mm]', 'daily averaged vapor pressure deficit [hpa]', 'daily saturated vapor pressure using daily averaged air temperature [hpa]', 'daily averaged relative humidity calculated [%]']
# varname = ["Air Temperature", "Net Radiation", "Gamma", "Slope of Saturation Vapor Pressure Curve", "Vapor Pressure Deficit", "Wind", "Latent Heat Constant", "Saturated Vapor Pressure", "Humidity", 'Outgoing Shortwave Radiation', 'Outgoing Longwave Radiation']

variables = ['tavg', 'sun', 'wind', 'ea', 'lc', 'ts', 'vpd', 'rain']
vars_penpan = ['tavg', 'tmax', 'tmin', 'p', 'ea', 'wind', 'tc', 'lat', 'elev'] #'sun'
basinlongs=['Songhuajiang','Liaohe','Northwestern','Haihe','Yellow', 'Yangtze','Huaihe','Southeastern','Southwestern','Pearl']
geoinfo = load('%s/station_geoinfo' %workspace)
station_number = load('%s/basin_station_number' %workspace)
# Dimensions
styr = 1961
edyr = 2001
stdy = datetime.datetime(styr, 1, 1)
eddy = datetime.datetime(edyr, 12, 31)
dates = pd.date_range(stdy, eddy, freq='D')
tstep = len(dates)
dyears = dates.year
dmonths = dates.month
doys = vstack([dates[i].timetuple().tm_yday for i in xrange(0, tstep)])

## quality check using data availability as criteria
station_flag = pickle.load(open('station_flag','rb'))
station_qc = [np.where(station_flag[ibasin][:, 0]==0)[0] for ibasin in xrange(0, 10)]
station_pan_flag = pickle.load(open('station_pan_flag','rb'))
station_pan_qc = [np.where(station_pan_flag[ibasin][:, 0]==0)[0] for ibasin in xrange(0, 10)]
good_stations = [intersect1d(station_qc[i], station_pan_qc[i]) for i in xrange(0, 10)]
##########################################################################
# for spectral coherece analysis
##########################################################################
nf = 513
nbasin = 10
nvar = 8
sampling_frequency = 1/(24.0 * 3600.0)  # unit: per day

def Coherence_Frequency():

	data = scipy.io.loadmat('%s/1_AP.mat' %(datadir))
	input = data[variables[0]][0, 0].flatten()
	pan = data['pan'][0, 0].flatten()
	freq = FFT.Coherence(input, pan, sampling_frequency, 'linear')[0]

	return freq

def Coherence_Station():

	for ibasin in xrange(1, 2): #11):
		for istation in xrange(0, 1):  #station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin))
			plt.plot(data['tavg'][0,0]-(data['tmax'][0,0]+data['tmin'][0,0])/2)
			plt.show()
			exit()
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [data[v][0, 0].flatten() for v in variables]
			pan = data['pan'][0, 0].flatten()
			# result = FFT.Power_Spectrum(pan, sampling_frequency, 'linear')
			coh = [FFT.Coherence(var, pan, sampling_frequency, 'linear')[1] for var in input]
			freq = FFT.Coherence(input[0], pan, sampling_frequency, 'linear')[0]
			Plotting.CoherenceStationPlot(coh, sampling_frequency, freq)

	return
## average the coherence spectrum
def Coherence_Station2Basin():

	for ibasin in xrange(0, 10):
		cohere_basin = []
		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			## Quality check
			input = [np.isnan(nanmean(data[v][0, istation])) for v in variables]
			flag = [1 for i in input if i == True]
			if len(flag) > 0:
				print "Ignore %s station %s!" % (basinlongs[ibasin], istation)
				continue

			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [data[v][0, istation].flatten() for v in variables]
			pan = data['pan'][0, istation].flatten()
			# Compute the coherence
			cohere_basin.append(vstack([FFT.Coherence(v, pan, sampling_frequency, 'linear')[1] for v in input]).reshape(1, nvar, nf))

		# store basin average
		cohere_basin = vstack(cohere_basin)
		cohere_basin.dump('%s/coherence_allvar_%s' %(workspace, basinlongs[ibasin]))

	return

def PSD_Station2Basin():

	for ibasin in xrange(0, 10):
		psd_basin = []
		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			## Quality check
			input = [np.isnan(nanmean(data[v][0, istation])) for v in variables]
			flag = [1 for i in input if i == True]
			if len(flag) > 0:
				print "Ignore %s station %s!" % (basinlongs[ibasin], istation)
				continue

			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [data[v][0, istation].flatten() for v in variables]
			# pan = data['pan'][0, istation].flatten()
			# Compute the coherence
			psd_basin.append(vstack([FFT.Power_Spectrum(v, sampling_frequency, 'linear')[1] for v in input]).reshape(1, nvar, nf))

		# store basin average
		psd_basin = vstack(psd_basin)
		psd_basin.dump('%s/psd_allvar_%s' %(workspace, basinlongs[ibasin]))

	return

def PlotCoherence():

	fig = plt.figure(figsize=(24, 18))
	cohere = []
	freq = Coherence_Frequency()
	for ibasin in xrange(0, 10):
		cohere_basin = load('%s/coherence_allvar_%s' %(workspace, basinlongs[ibasin]))
		cohere.append(mean(cohere_basin, axis=0).reshape(1, nvar, nf))
		## DRAW FIGURE---------------
	 	ax = fig.add_subplot(3, 4, ibasin+1)
		Plotting.CoherenceBasinPlot(ax, mean(cohere_basin, axis=0), sampling_frequency, freq, basinlongs[ibasin])

	# for national average
	ax = fig.add_subplot(3, 4, 11)
	Plotting.CoherenceBasinPlot(ax, mean(vstack(cohere), axis=0), sampling_frequency, freq, 'Average')
	ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center left', fontsize=21)
	fig.tight_layout()
	plt.show()

def PlotPSD():

	fig = plt.figure(figsize=(24, 18))
	psd = []
	freq = Coherence_Frequency()
	for ibasin in xrange(0, 10):
		psd_basin = load('%s/psd_allvar_%s' %(workspace, basinlongs[ibasin]))
		psd.append(mean(psd_basin, axis=0).reshape(1, nvar, nf))
		## DRAW FIGURE---------------
	 	ax = fig.add_subplot(3, 4, ibasin+1)
		Plotting.PSDBasinPlot(ax, mean(psd_basin, axis=0), sampling_frequency, freq, basinlongs[ibasin])

	# for national average
	ax = fig.add_subplot(3, 4, 11)
	Plotting.PSDBasinPlot(ax, mean(vstack(psd), axis=0), sampling_frequency, freq, 'Average')
	ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center left', fontsize=21)
	fig.tight_layout()
	plt.show()


##########################################################################
# trend attribution analysis for pan evporation and penpan modelled PE
##########################################################################

###==================================
# Step1. Select good station: no missing years and missing months
def Quality_Check_station():
	"calculate the total number of missing records for each station"
	flag_basins = []
	for ibasin in xrange(0, 10):
		flag = zeros((station_number[ibasin], 2))
		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			# for a station to count it must have more than 80% of coverage for each year and for all variables for Epan calculation
			for v in vars_penpan[:-2]:
				ts = data[v][0, istation][0:tstep]

				# Case1: If the whole series is 0, then throw away
				if np.nansum(ts) == 0.0:
					flag[istation, 0] = 1
					flag[istation, 1] = 1
					break
				# # Case2: For every month, if the missing value is greater than the threshold 20 % then mask down, stricter than year
				for year in xrange(styr, edyr+1):
					for month in xrange(1, 13):
						npoints_year = 31.0*(1-80.0/100.0)

						count = np.isnan(ts[(dyears==year)&(dmonths==month)]).sum()

						if count > int(npoints_year):
							flag[istation, 0] = 1
							flag[istation, 1] = 2
							break

		flag_basins.append(flag)

	return flag_basins

# station_flag = Quality_Check_station()
# pickle.dump(station_flag, open('station_flag','wb'))

def Quality_Check_station_pan():
	"calculate the total number of missing records for each station"
	flag_basins = []
	for ibasin in xrange(0, 10):
		flag = zeros((station_number[ibasin], 2))
		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			# for a station to count it must have more than 80% of coverage for each year and for all variables for Epan calculation
			ts = data['pan'][0, istation][0:tstep]
			# Case1: If the whole series is 0, then throw away
			if np.nansum(ts) == 0.0:
				flag[istation, 0] = 1
				flag[istation, 1] = 1
				break
			# # Case2: For every month, if the missing value is greater than the threshold 20 % then mask down, stricter than year
			for year in xrange(styr, edyr+1):
				for month in xrange(1, 13):
					npoints_year = 31.0*(1-80.0/100.0)
					count = np.isnan(ts[(dyears==year)&(dmonths==month)]).sum()
					if count > int(npoints_year):
						flag[istation, 0] = 1
						flag[istation, 1] = 2
						break

		flag_basins.append(flag)

	return flag_basins
# station_pan_flag = Quality_Check_station_pan()
# pickle.dump(station_pan_flag, open('station_pan_flag','wb'))

## use the R square between pan evaporation and Epan to filter out some stations further
def Compare_pan_vs_Epenpan(INPUT, pan):

	test = Data(INPUT, 'cloud').penpan
	res = pd.Series(pan, index=dates).resample('M', how='sum')
	res2 = pd.Series(test, index=dates).resample('M', how='sum')
	plt.scatter(res, res2)
	mk = np.isnan(res)+np.isnan(res2)
	# calculate R square
	R = scipy.stats.mstats.pearsonr(res[mk==0], res2[mk==0])[0]
	plt.title('%s' %(R**2))
	plt.xlim([0, 500])
	plt.ylim([0, 500])
	plt.xlabel('observed')
	plt.ylabel('penpan')
	plt.show()
	exit()

	return

###==================================
# Step2. Gapfill the missing values in the selected stations
def Gapfill(daily):
	"return the array, not pandas object"
	ts = pd.Series(daily, index=dates).fillna(method='pad').values
	return ts

## Step4. Calculate Epan trend attribution
def generate_annual_trend_attribution_components(INPUT):

	"""prepare all the trend attribution components"""
	# Net Radiation Rn:  (delta_ave, gamma_ave, trend_Rnet)
	# Wind speed U:  (delta_ave, gamma_ave, vpd_ave, lhv_ave, trend_Wind)
	# Humidity Ea: (gamma_ave, wind_ave, delta_ave, lhv_ave, trend_Eact)
	# Temperature Ta: (Rn_ave, gamma_ave, delta_ave, trend_Delta, vpd_ave, wind_ave, lhv_ave, trend_Esat)
	res = Data(INPUT, 'cloud')
	Tair_daily = res.Tair
	Rn_daily = res.Rn_pan
	PDRnl_daily = res.PDLWnet_PDTair
	Gamma_daily = res.gamma
	Delta_daily = res.DELTA
	Vpd_daily = res.vpd
	Wind_daily = res.Wind
	Lhv_daily = res.lv
	Es_daily = res.estar
	Ea_daily = res.e
	del res

	vars = [Tair_daily, Rn_daily, PDRnl_daily, Gamma_daily, Delta_daily, Vpd_daily, Wind_daily, Lhv_daily, Es_daily, Ea_daily]


	# aggregate the time series into annual time series for trend analysis
	def daily2annual(daily):
		ts = pd.Series(daily, index=dates).fillna(method='pad')
		return ts.resample('A', how='mean')
	forcing_ann = [daily2annual(var) for var in vars]

	return forcing_ann

def Calculate_forcing_ave_trend4trend(forcing_ann):

	"""calculate ave and trend for the trend attribution components"""
	# Net Radiation Rn:  (delta_ave, gamma_ave, trend_Rnet)
	# Wind speed U:  (delta_ave, gamma_ave, vpd_ave, lhv_ave, trend_Wind)
	# Humidity Ea: (gamma_ave, wind_ave, delta_ave, lhv_ave, trend_Eact)
	# Temperature Ta: (Tair_ave, Rn_ave, gamma_ave, delta_ave, trend_Delta, trend_Tair, vpd_ave, wind_ave, lhv_ave, trend_Esat)
	# vars = [Tair_daily, Rn_daily, PDRnl_daily, Gamma_daily, Delta_daily, Vpd_daily, Wind_daily, Lhv_daily, Es_daily, Ea_daily]
	avenames = ['delta_ave', 'gamma_ave', 'vpd_ave', 'lhv_ave', 'Rn_ave', 'wind_ave', 'PDRnl_ave']
	avelist = [4, 3, 5, 7, 1, 6, 2]
	trendnames = ['trend_Rn', 'trend_Wind', 'trend_Eact', 'trend_Delta', 'trend_Tair', 'trend_Esat']
	trendlist = [1, 6, 9, 4, 0, 8]
	forave = {avenames[i]: forcing_ann[avelist[i]].mean() for i in xrange(0, len(avelist))}
	fortrend = {trendnames[i]: TrendAnalysis(forcing_ann[trendlist[i]]).ThSenTrend() for i in xrange(0, len(trendlist))}

	return forave, fortrend

def Calculate_trend_attribution(forcing_ann):

	"""Calculate the PE trend attribution component"""
	forave, fortrend = Calculate_forcing_ave_trend4trend(forcing_ann)
	attribution = TrendAttribution(forave, fortrend).Calculate_Attribution()

	return attribution

##====================================================================
# Run the step functions for spectral analysis
##====================================================================
# Coherence_Station()
# Coherence_Station2Basin()
# PlotCoherence()
# PSD_Station2Basin()
# PlotPSD()

##====================================================================
# Run the step functions for trend attribution analysis
##====================================================================
def Calculate_Pan():
	lons = []
	lats = []
	Pan_obs = []
	Pan_obs_gapfill = []
	Pan_mod = []
	for ibasin in xrange(0, 10):
		for istation in good_stations[ibasin]:
			print ibasin, istation
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
			pan_obs = data['pan'][0, istation][0:tstep].flatten()
			pan_obs_gapfill = Gapfill(data['pan'][0, istation][0:tstep].flatten())
			## Prepare for the input data for Epan calculation
			INPUT = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
			INPUT['doy'] = doys.flatten()
			INPUT['lat'] = geoinfo[index, 1]
			INPUT['elev'] = geoinfo[index, 3]

			## Step1. Calculate Epan and evaluation
			pan_mod = Data(INPUT, 'cloud').penpan
			## Step2. Calculate Epan Trend attribution
			# forcing_ann = generate_annual_trend_attribution_components(INPUT)
			# att = Calculate_trend_attribution(forcing_ann)

			# Collect the data and geoinfo for map showing
			lons.append(geoinfo[index, 2])
			lats.append(geoinfo[index, 1])
			Pan_obs.append(pan_obs)
			Pan_obs_gapfill.append(pan_obs_gapfill)
			Pan_mod.append(pan_mod)
			# R square
			#
			# print att
	lons = vstack(lons)
	lats = vstack(lats)
	Pan_obs = vstack(Pan_obs)
	Pan_obs_gapfill = vstack(Pan_obs_gapfill)
	Pan_mod = vstack(Pan_mod)
	lons.dump('%s/lons_good_stations' %workspace)
	lats.dump('%s/lats_good_stations' %workspace)
	Pan_obs.dump('%s/pan_obs_ts_good_stations' %workspace)
	Pan_obs_gapfill.dump('%s/pan_obs_gapfill_good_stations' %workspace)
	Pan_mod.dump('%s/pan_mod_good_stations' %workspace)
	return

# Calculate_Pan()

## Plot the
# lons = load('%s/lons_good_stations' %workspace)
# lats = load('%s/lats_good_stations' %workspace)
# data = load('%s/pan_obs_ts_good_stations' %workspace)
# data = nanmean(data, axis=1)
# print data.shape
# lons = []
# lats = []
# for ibasin in xrange(0, 10):
# 	for istation in xrange(0, station_number[ibasin]):
# 		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
# 		index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
# 		lons.append(geoinfo[index, 2])
# 		lats.append(geoinfo[index, 1])
# lons = vstack(lons)
# lats = vstack(lats)
# lons.dump('%s/lons_all_stations' %workspace)
# lats.dump('%s/lats_all_stations' %workspace)
lons = load('%s/lons_all_stations' %workspace)
lats = load('%s/lats_all_stations' %workspace)
Plotting.Mapshow(ones(lons.shape), lons, lats, 0, 10, plt.cm.jet, None, None, figdir, None)

exit()
# Plotting.Mapshow(data, lons, lats, 0, 10, plt.cm.jet, None, None, figdir, None)

