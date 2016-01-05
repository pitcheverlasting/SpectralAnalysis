__author__ = 'lpeng'

from pylab import *
import pickle
import scipy.io
import scipy.stats
from scipy.optimize import minimize
import statsmodels.tools.eval_measures as evaluate
import pandas as pd
import IO, FFT, Plotting
from PET_Library import Data
from PETREND import TrendAnalysis, TrendAttribution
# #=============================path================================
datadir = '/home/water5/lpeng/script/PROJECT/pan_evaporation_china/original'
workspace = '/home/water5/lpeng/script/PROJECT/pan_evaporation_china/workspace/201601'
figdir = '/home/water5/lpeng/Figure/pan_spectral/sunhours'
##=============================variable============================
vars = ['time', 'p', 'tavg', 'tmin', 'tmax', 'ea', 'rh', 'tc', 'lc', 'wind', 'ts', 'sun', 'rain', 'pan', 'vpd',  'estavg', 'rh_test']  # time 1,2,3 year, month, day # Note that original data is multiplied by 10 from the normal unit. This process data devide 10 so that it becomes the original unit
varlongs = ['time', 'daily averaged pressure [hpa]', 'daily averaged air temperature [degree]', 'daily minimum air temperature [degree]', 'daily maximum air temperature [degree]', 'daily averaged vapor pressure [hpa]', 'daily averaged relative humidity [%]', 'daily averaged total cloud cover', 'daily averaged low cloud cover', 'daily averaged wind speed [m/s]', 'daily averaged surface temperature(at 0cm level) [degree]', 'daily averaged sunshine hour [hr]', 'daily averaged rainfall [mm]', 'daily averaged pan evaporation [mm]', 'daily averaged vapor pressure deficit [hpa]', 'daily saturated vapor pressure using daily averaged air temperature [hpa]', 'daily averaged relative humidity calculated [%]']

variables = ['tavg', 'sun', 'wind', 'ea', 'tc', 'vpd']
vars_penpan = ['tavg', 'tmax', 'tmin', 'p', 'ea', 'wind', 'sun', 'lat', 'elev']  # 'tc'
basinlongs=['Songhuajiang', 'Liaohe', 'Northwestern', 'Haihe', 'Yellow', 'Yangtze', 'Huaihe', 'Southeastern', 'Southwestern', 'Pearl']
geoinfo = load('%s/station_geoinfo' %workspace)
station_number = load('%s/basin_station_number' %workspace)

## Time
styr = 1961
edyr = 2001
stdy = datetime.datetime(styr, 1, 1)
eddy = datetime.datetime(edyr, 12, 31)
dates = pd.date_range(stdy, eddy, freq='D')
tstep = len(dates)
dyears = dates.year
dmonths = dates.month
doys = vstack([dates[i].timetuple().tm_yday for i in xrange(0, tstep)])  # julian day

## quality check using data availability as criteria
station_flag = pickle.load(open('station_sunhours_80_flag','rb'))
station_qc = [np.where(station_flag[ibasin][:, 0]==0)[0] for ibasin in xrange(0, 10)]
station_pan_flag = pickle.load(open('station_pan_80_flag','rb'))
station_pan_qc = [np.where(station_pan_flag[ibasin][:, 0]==0)[0] for ibasin in xrange(0, 10)]
good_stations = [intersect1d(station_qc[i], station_pan_qc[i]) for i in xrange(0, 10)]

##########################################################################
# construct database
##########################################################################
# Step1. Select good station: no missing years and missing months 80% coverage
# Step1.1 Select good met station
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
				npoints_mon = 31.0*(1-80.0/100.0)
				for year in xrange(styr, edyr+1):
					for month in xrange(1, 13):
						count = np.isnan(ts[(dyears==year)&(dmonths==month)]).sum()

						if count > int(npoints_mon):
							flag[istation, 0] = 1
							flag[istation, 1] = 2
							break

		flag_basins.append(flag)

	return flag_basins

# station_flag = Quality_Check_station()
# pickle.dump(station_flag, open('station_sunhours_80_flag', 'wb'))
# exit()

# Step1.2 Select good pan station
def Quality_Check_station_pan(percentage):
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
					npoints_year = 31.0*(1-percentage/100.0)
					count = np.isnan(ts[(dyears==year)&(dmonths==month)]).sum()
					if count > int(npoints_year):
						flag[istation, 0] = 1
						flag[istation, 1] = 2
						break

		flag_basins.append(flag)

	pickle.dump(station_pan_flag, open('station_pan_%2d_flag','wb') % percentage)

	return flag_basins

# station_pan_flag = Quality_Check_station_pan(60)
# exit()

station_list_names = ['_all_stations', '_good_stations', '_for_stations', '_pan_stations']
station_lists = [station_number, good_stations, station_qc, station_pan_qc]
fignames = ['all', 'good', 'forcing', 'pan']
cols = ['grey', 'springgreen', 'aquamarine', 'lightcoral']

# Step1.3 Plot good stations
def Get_Selected_Station_Locations(j):
	lons = []
	lats = []
	for ibasin in xrange(0, 10):
		for istation in station_lists[j][ibasin]:
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
			lons.append(geoinfo[index, 2])
			lats.append(geoinfo[index, 1])
	lons = vstack(lons)
	lats = vstack(lats)
	lons.dump('%s/lons%s' %(workspace, station_list_names[j]))
	lats.dump('%s/lats%s' %(workspace, station_list_names[j]))

	return

def Plot_Selected_Station_Locations(j):
	lons = load('%s/lons%s' %(workspace, station_list_names[j]))
	lats = load('%s/lats%s' %(workspace, station_list_names[j]))
	Plotting.Mapshow(cols[j], lons, lats, 45, 0.5, None, None, None, 'Qualified pan evaporation met stations (%s)' %len(lons), 'Met Stations', '', figdir, '%s_stations_scattermap.png' %fignames[j])

	return

# Get_Selected_Station_Locations(1)
# Plot_Selected_Station_Locations(1)
# exit()

###====================================================================
# Step1.4 Gapfill the missing values in the selected stations
def Gapfill(daily):
	"return the array, not pandas object"
	ts = pd.Series(daily, index=dates).fillna(method='pad').values
	return ts

def daily2annual(daily):
	ts = pd.Series(daily, index=dates).resample('A', how='mean').values
	return ts

def daily2monthly(daily):
	ts = pd.Series(daily, index=dates).resample('M', how='mean').values
	return ts
###====================================================================

##########################################################################
# ## parameterization
##########################################################################
def KGE(model, obs):

	rho = scipy.stats.pearsonr(model, obs)
	mean_ratio = np.mean(model) / np.mean(obs)
	std_ratio = np.std(model) / np.std(obs)

	return 1 - ((rho[0] - 1) ** 2 + (mean_ratio - 1) ** 2 + (std_ratio - 1) ** 2) ** 0.5

def Calibrate_fq():
	def Function(params):
		a, b = params
		res = []
		for ibasin in xrange(0, 1): #10):
			for istation in good_stations[ibasin]:
				# print ibasin, istation
				data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
				index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
				# pan_obs = data['pan'][0, istation][0:tstep].flatten()
				pan_obs_gapfill = Gapfill(data['pan'][0, istation][0:tstep].flatten())
				## Prepare for the input data for Epan calculation
				INPUT = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
				INPUT['doy'] = doys.flatten()
				INPUT['lat'] = geoinfo[index, 1]
				INPUT['elev'] = geoinfo[index, 3]
				pan_mod = Data(INPUT, 'cloud').Penpan_u2(a, b)
				res.append(evaluate.rmse(daily2monthly(pan_mod), daily2monthly(pan_obs_gapfill)))

		return vstack(res).mean()
			# exit()
			# A.append(res['x'][0])
			# B.append(res['x'][1])
	bnd = ((0, None), (0, None))
	test = minimize(Function, array([1.313, 1.381]), bounds=bnd, method='SLSQP')
	print test


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
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [data[v][0, 0].flatten() for v in variables]
			pan = data['pan'][0, 0].flatten()
			# result = FFT.Power_Spectrum(pan, sampling_frequency, 'linear')
			coh = [FFT.Coherence(var, pan*0.9+0.1, sampling_frequency, 'linear')[1] for var in input]
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

##====================================================================
# Run the step functions for spectral analysis
##====================================================================
# Coherence_Station()
# Coherence_Station2Basin()
# PlotCoherence()
# PSD_Station2Basin()
# PlotPSD()


##########################################################################
# Calculate penpan modelled PE
##########################################################################

def Calculate_Epenpan_daily():
	Pan_obs_gapfill = []
	Pan_mod = []
	for ibasin in xrange(0, 10):
		for istation in good_stations[ibasin]:
			print ibasin, istation
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
			pan_obs_gapfill = Gapfill(data['pan'][0, istation][0:tstep].flatten())
			## Prepare for the input data for Epan calculation
			INPUT = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
			INPUT['doy'] = doys.flatten()
			INPUT['lat'] = geoinfo[index, 1]
			INPUT['elev'] = geoinfo[index, 3]

			## Calculate Epan and evaluation
			pan_mod = Data(INPUT, 'sunhours').penpan
			# Collect the data and geoinfo for map showing
			Pan_obs_gapfill.append(pan_obs_gapfill)
			Pan_mod.append(pan_mod)
	Pan_obs_gapfill = vstack(Pan_obs_gapfill)
	Pan_mod = vstack(Pan_mod)
	Pan_obs_gapfill.dump('%s/pan_obs_gapfill_good_stations' %workspace)
	Pan_mod.dump('%s/pan_mod_good_stations' %workspace)

	return

# Calculate_Epenpan_daily()
# exit()

def Calculate_radiation_daily():
	Rs = []
	Rnl = []
	Rn = []
	for ibasin in xrange(0, 10):
		for istation in good_stations[ibasin]:
			print ibasin, istation
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
			## Prepare for the input data for Epan calculation
			INPUT = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
			INPUT['doy'] = doys.flatten()
			INPUT['lat'] = geoinfo[index, 1]
			INPUT['elev'] = geoinfo[index, 3]

			## Calculate
			result = Data(INPUT, 'sunhours')
			Rs.append(result.SWin)
			Rnl.append(result.LWnet)
			Rn.append(result.Rn_pan)
	Rs = vstack(Rs)
	Rnl = vstack(Rnl)
	Rn = vstack(Rn)
	Rs.dump('%s/Rs_0.25_0.5_good_stations' %(workspace))
	Rnl.dump('%s/Rnl_FAO_0.25_0.5_good_stations' %(workspace))
	Rn.dump('%s/Rn_0.25_0.5_good_stations' %(workspace))

	return
# Calculate_radiation_daily()
# exit()

def Coherence_obs_radiation():

	## Prepare the data for plotting
	# pan_obs_gapfill = load('%s/pan_obs_gapfill_qc_stations' %workspace)
	r2_60_stations = pickle.load(open('r2_60_stations','rb'))
	r2_60 = hstack(r2_60_stations)
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)[r2_60,:]
	Rs = load('%s/Rs_0.25_0.5_good_stations' %(workspace))[r2_60,:]
	Rnl = load('%s/Rnl_FAO_0.25_0.5_good_stations' %(workspace))[r2_60,:]
	Rn = load('%s/Rn_0.25_0.5_good_stations' %(workspace))[r2_60,:]

	relatetable = load('%s/relatetable' %workspace)
	ist = 0
	for ibasin in xrange(0, 10):
		cohere_obs_basin = []
		for istation in relatetable[r2_60_stations[ibasin], 1]: #station_pan_qc[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			panobs = pan_obs_gapfill[ist, :].flatten()
			rs = Rs[ist, :].flatten()
			rnl = Rnl[ist, :].flatten()
			rn = Rn[ist, :].flatten()
			rad = [rs, rnl, rn]
			# Compute the coherence
			cohere_obs_basin.append(vstack([FFT.Coherence(r, panobs, sampling_frequency, 'linear')[1] for r in rad]).reshape(1, 3, nf))
			ist = ist + 1

		cohere_obs_basin = vstack(cohere_obs_basin)
		cohere_obs_basin.dump('%s/coherence_obs_radiation_0.25_0.5_r2_60_station_%s' %(workspace, basinlongs[ibasin]))

# Coherence_obs_radiation()
# exit()


def Coherence_obs_qc():

	## Prepare the data for plotting
	# pan_obs_gapfill = load('%s/pan_obs_gapfill_qc_stations' %workspace)

	r2_60_stations = pickle.load(open('r2_60_stations','rb'))
	r2_60 = hstack(r2_60_stations)
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)[r2_60,:]
	relatetable = load('%s/relatetable' %workspace)
	ist = 0
	for ibasin in xrange(0, 10):
		cohere_obs_basin = []
		for istation in relatetable[r2_60_stations[ibasin], 1]: #station_pan_qc[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			input = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in variables]
			panobs = pan_obs_gapfill[ist, :].flatten()
			# Compute the coherence
			cohere_obs_basin.append(vstack([FFT.Coherence(v, panobs, sampling_frequency, 'linear')[1] for v in input]).reshape(1, nvar, nf))
			ist = ist + 1

		# store basin average
		cohere_obs_basin = vstack(cohere_obs_basin)
		# cohere_obs_basin.dump('%s/coherence_obs_var_qc_station_%s' %(workspace, basinlongs[ibasin]))
		cohere_obs_basin.dump('%s/coherence_obs_var_r2_60_station_%s' %(workspace, basinlongs[ibasin]))

	return
# Coherence_obs_qc()
# exit()


def Coherence_mod_r260():

	## Prepare the data for plotting
	r2_60_stations = pickle.load(open('r2_60_stations_u56','rb'))
	r2_60 = hstack(r2_60_stations)
	pan_mod = load('%s/pan_mod_good_stations_u56' %workspace)[r2_60,:]
	relatetable = load('%s/relatetable' %workspace)
	stats_basin = load('%s/calibration_regressionfit_basin_u56' %workspace)

	ist = 0
	for ibasin in xrange(0, 10):
		cohere_mod_basin = []
		for istation in relatetable[r2_60_stations[ibasin], 1]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			input = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in variables]
			panmod = pan_mod[ist, :].flatten() * stats_basin[ibasin, 0] #+ stats_basin[ibasin, 1]
			# Compute the coherence
			cohere_mod_basin.append(vstack([FFT.Coherence(v, panmod, sampling_frequency, 'linear')[1] for v in input]).reshape(1, nvar, nf))
			ist = ist + 1

		# store basin average
		cohere_mod_basin = vstack(cohere_mod_basin)
		cohere_mod_basin.dump('%s/coherence_mod_var_goodstation_r2_60_u56_%s' %(workspace, basinlongs[ibasin]))

	return
# Coherence_mod_r260()
# exit()

def Coherence_obs_mod_r260():

	## Prepare the data for plotting
	r2_60_stations = pickle.load(open('r2_60_stations_u56','rb'))
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	pan_mod = load('%s/pan_mod_good_stations' %workspace)
	stats_basin = load('%s/calibration_regressionfit_basin_u56' %workspace)
	for ibasin in xrange(0, 10):
		cohere_obs_mod_basin = []
		for istation in good_stations[ibasin]: #r2_60_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			panobs = pan_obs_gapfill[istation, :].flatten()
			# panmod = pan_mod[istation, :].flatten()
			panmod = pan_mod[istation, :].flatten() * stats_basin[ibasin, 0]
			# Compute the coherence
			# data = FFT.Coherence(panobs, panobs, sampling_frequency, 'linear')[1]
			# Plotting.CoherenceStationPlot(data, sampling_frequency, Coherence_Frequency())
			# exit()
			# cohere_obs_mod_basin.append(FFT.Coherence(panobs, panmod, sampling_frequency, 'linear')[1])

		# store basin average
		cohere_obs_mod_basin = vstack(cohere_obs_mod_basin)
		cohere_obs_mod_basin.dump('%s/coherence_obs_mod_var_goodstation_r2_60_u56_%s' %(workspace, basinlongs[ibasin]))

	return
# Coherence_obs_mod_r260()
# exit()

def Plot_Coherence_Basin():

	fig = plt.figure(figsize=(24, 18))
	cohere = []
	freq = Coherence_Frequency()
	for ibasin in xrange(0, 10):
		# cohere_basin = load('%s/coherence_obs_var_qc_station_%s' %(workspace, basinlongs[ibasin]))
		# cohere_basin = load('%s/coherence_mod_var_goodstation_%s' %(workspace, basinlongs[ibasin]))
		# cohere_basin = load('%s/coherence_obs_var_r2_60_station_%s' %(workspace, basinlongs[ibasin]))
		# cohere_basin = load('%s/coherence_mod_var_goodstation_r2_60_u56_%s' %(workspace, basinlongs[ibasin]))
		cohere_basin = load('%s/coherence_obs_radiation_%s_r2_60_station_%s' %(workspace, 'cloud', basinlongs[ibasin]))
		cohere.append(cohere_basin)

		## DRAW FIGURE---------------
	 	ax = fig.add_subplot(3, 4, ibasin+1)
		Plotting.CoherenceBasinPlot(ax, mean(cohere_basin, axis=0), sampling_frequency, freq, basinlongs[ibasin])

	# for national average
	ax = fig.add_subplot(3, 4, 11)
	Plotting.CoherenceBasinPlot(ax, mean(vstack(cohere), axis=0), sampling_frequency, freq, 'Average')
	ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center left', fontsize=21)
	fig.tight_layout()
	plt.show()
	return
Plot_Coherence_Basin()
exit()

def Plot_Coherence_obs_mod():

	fig = plt.figure(figsize=(24, 18))
	cohere = []
	freq = Coherence_Frequency()
	for ibasin in xrange(0, 10):
		cohere_basin = load('%s/coherence_obs_mod_var_goodstation_r2_60_u56_%s' %(workspace, basinlongs[ibasin]))
		# cohere_basin = load('%s/coherence_obs_mod_var_goodstation_%s' %(workspace, basinlongs[ibasin]))
		cohere.append(cohere_basin)

		## DRAW FIGURE---------------
	 	ax = fig.add_subplot(3, 4, ibasin+1)
		# Plotting.CoherenceBasinPlot(ax, mean(cohere_basin, axis=0), sampling_frequency, freq, basinlongs[ibasin])

		Plotting.Coherence_obs_mod_Plot(ax, mean(cohere_basin, axis=0), sampling_frequency, freq, basinlongs[ibasin])

	# for national average
	ax = fig.add_subplot(3, 4, 11)
	# Plotting.CoherenceBasinPlot(ax, mean(vstack(cohere), axis=0), sampling_frequency, freq, 'Average')
	Plotting.Coherence_obs_mod_Plot(ax, mean(vstack(cohere), axis=0), sampling_frequency, freq, 'Average')
	fig.tight_layout()
	plt.show()
# Plot_Coherence_obs_mod()
# exit()
##########################################################################
# trend attribution analysis for pan evporation and penpan modelled PE
##########################################################################

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

def Epenpan_trend_attribution():

	att_stations = []
	r2_60_stations = pickle.load(open('r2_60_stations_u56','rb'))
	relatetable = load('%s/relatetable' %workspace)

	for ibasin in xrange(0, 10):
		for istation in relatetable[r2_60_stations[ibasin], 1]:
			print ibasin, istation
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))

			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]

			## Prepare for the input data for Epan calculation
			INPUT = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
			INPUT['doy'] = doys.flatten()
			INPUT['lat'] = geoinfo[index, 1]
			INPUT['elev'] = geoinfo[index, 3]

			## Calculate Epan Trend attribution
			forcing_ann = generate_annual_trend_attribution_components(INPUT)
			att = Calculate_trend_attribution(forcing_ann)
			# Collect the data and geoinfo for map showing
			att_stations.append(att)

	att_stations = vstack(att_stations)
	att_stations.dump('%s/att_stations' %workspace)

	return

# Epenpan_trend_attribution()
# exit()


##====================================================================
# Run the step functions for trend attribution analysis
##====================================================================
def Compare_pan_vs_Epenpan():
	## Prepare the data for plotting
	# pan_obs = load('%s/pan_obs_ts_good_stations' %workspace)
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	pan_mod = load('%s/pan_mod_good_stations_u56' %workspace)

	stats_basin = load('%s/calibration_regressionfit_basin_u56' %workspace)
	r2_60_stations = pickle.load(open('r2_60_stations_u56','rb'))
	r2_60 = hstack(r2_60_stations)

	# station average
	# pan_obs_gapfill = mean(vstack([daily2monthly(pan_obs_gapfill[istation, :]) for istation in xrange(0, 236)]), axis=0)
	# pan_mod = mean(vstack([daily2monthly(pan_mod[istation, :]) for istation in xrange(0, 236)]), axis=0)
	# all stations all monthly time step
	# panobs = pan_obs_gapfill[r2_60,:]
	# panobs = vstack([daily2monthly(panobs[istation, :]) for istation in xrange(0, len(r2_60))]).flatten()
	# panmod = []
	# for ibasin in xrange(0, 10):
	# 	for istation in r2_60_stations[ibasin]:
	# 		panmod.append(daily2monthly(pan_mod[istation, :]) * stats_basin[ibasin, 0]) #  + stats_basin[ibasin, 1]
	# panmod = vstack(panmod).flatten()
	# stats = scipy.stats.mstats.linregress(panobs, panmod)[:]
	# rmse = np.sqrt(mean((panmod-panobs)**2))

	## annual trend
	# all stations
	# panobs = vstack([TrendAnalysis(daily2annual(pan_obs_gapfill[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, 236)]).flatten()
	# panmod = vstack([TrendAnalysis(daily2annual(pan_mod[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, 236)]).flatten()
	# only high agreement stations

	panobs = pan_obs_gapfill[r2_60,:]
	panmod = []
	for ibasin in xrange(0, 10):
		for istation in r2_60_stations[ibasin]:
			panmod.append(pan_mod[istation, :].flatten() * stats_basin[ibasin, 0]) # + stats_basin[ibasin, 1])
	panmod = vstack(panmod)
	panobs = vstack([TrendAnalysis(daily2annual(panobs[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, len(r2_60))]).flatten()
	panmod = vstack([TrendAnalysis(daily2annual(panmod[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, len(r2_60))]).flatten()
	#
	# mk_same_sign = np.where((panobs<0)&(panmod<0))[0]
	# plt.scatter(panobs, panmod, s=100, marker='+', c='lightCoral', linewidths=2, alpha=0.2) #
	fig = plt.figure(figsize=(8,8))
	plt.scatter(panobs, panmod, s=120, c='green', alpha=0.7) #

	# calculate R square
	R = scipy.stats.mstats.pearsonr(panobs, panmod)[0]
	min = -20
	max = 20
	plt.xlim([min, max])
	plt.ylim([min, max])
	plt.plot([min, max], [min, max], 'k--', linewidth=1.5)
	plt.xlabel('Observed (mm/yr/yr)', fontsize=21)
	plt.ylabel('Modelled (mm/yr/yr)', fontsize=21)
	plt.axes().set_aspect('equal')
	plt.axes().tick_params(axis='y', labelsize=19)
	plt.axes().tick_params(axis='x', labelsize=19)
	# plt.title('Monthly averaged pan evaporation (all stations)', fontsize=22, y=1.05)
	# plt.text(2, 18, r'$R^2$ = %0.3f' %(R**2), fontsize=20)
	plt.text(-18, 16, r'$R^2$ = %0.3f' %(R**2), fontsize=24)
	plt.title('Annual trend (all stations)', fontsize=16, y=1.03)
	plt.show()
	return

# Compare_pan_vs_Epenpan()
# exit()
def Calibrate_Epan_mod():
	# Calibrate the modelled pan evaporation for each basin
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	pan_mod = load('%s/pan_mod_good_stations_u56' %workspace)
	stats = [scipy.stats.mstats.linregress(pan_mod[istation, :], pan_obs_gapfill[istation, :])[:] for istation in xrange(0, 236)]
	stats = vstack(stats)
	# relatetable = []
	# for ibasin in xrange(0, 10):
	# 	for istation in good_stations[ibasin]:
	# 		relatetable.append((ibasin,istation))
	# relatetable = vstack(relatetable)
	# relatetable.dump('%s/relatetable' %workspace)
	relatetable = load('%s/relatetable' %workspace)
	stats_basin = []
	r2_60_stations = []
	for ibasin in xrange(0, 10):
		stations = np.where(relatetable[:, 0]==ibasin)[0]
		index = np.where(stats[stations, 2]**2>0.6)[0]
		r2_60_stations.append(stations[index])

		panobs = vstack([daily2monthly(pan_obs_gapfill[istation, :]) for istation in stations[index]]).flatten()
		panmod = vstack([daily2monthly(pan_mod[istation, :]) for istation in stations[index]]).flatten()
		stats_basin.append(scipy.stats.mstats.linregress(panmod, panobs)[:])
	stats_basin = vstack(stats_basin)
	stats_basin.dump('%s/calibration_regressionfit_basin_u56' %workspace)
	pickle.dump(r2_60_stations, open('r2_60_stations_u56', 'wb'))
	return

# Calibrate_Epan_mod()
# exit()
def Plot_Station_Data():

	stats_basin = load('%s/calibration_regressionfit_basin_u56' %workspace)
	r2_60_stations = pickle.load(open('r2_60_stations_u56','rb'))
	r2_60 = hstack(r2_60_stations)
	lons = load('%s/lons_good_stations' %workspace)[r2_60, :]
	lats = load('%s/lats_good_stations' %workspace)[r2_60, :]

	## Prepare the data for plotting
	# 1. mean
	# min = 0
	# max = 10
	# title = 'daily average'
	# unit = '(mm/d)'
	# figname = 'ave'
	# pan_obs_gapfill = mean(load('%s/pan_obs_gapfill_good_stations' %workspace), axis=1)
	# pan_mod = mean(load('%s/pan_mod_good_stations' %workspace), axis=1)

	# 2. inter annual variability
	# min = 0
	# max = 20
	# title = 'interannual coefficient of variability'
	# unit = '(%)'
	# figname = 'cv'
	# pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	# pan_obs_gapfill_ann = vstack([daily2annual(pan_obs_gapfill[istation, :]) for istation in xrange(0, len(lons))])
	# pan_obs_gapfill = scipy.stats.variation(pan_obs_gapfill_ann, axis=1) * 100.0
	#
	# pan_mod = load('%s/pan_mod_good_stations' %workspace)
	# pan_mod_ann = vstack([daily2annual(pan_mod[istation, :]) for istation in xrange(0, len(lons))])
	# pan_mod = scipy.stats.variation(pan_mod_ann, axis=1) * 100.0
	#
	# ##----------------------------------Figure Plotting-----------------------------------------
	# Plotting.Mapshow(pan_obs_gapfill, lons, lats, 50, min, max, plt.cm.jet, 'Observed pan evaporation %s' %title, None, unit, figdir, 'pan_obs_%s.png' % figname)
	# Plotting.Mapshow(pan_mod, lons, lats, 50, min, max, plt.cm.jet, 'Modelled pan evaporation %s' %title, None, unit, figdir, 'pan_mod_%s.png' % figname)

	# 3. annual trend
	# min = -25
	# max = 25
	# title = 'annual trend'
	# unit = '(mm/yr/yr)'
	# figname = 'trend'

	# pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)[r2_60,:]
	# panobs_slp = vstack([TrendAnalysis(daily2annual(pan_obs_gapfill[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, len(r2_60))])
	# panobs_p = vstack([TrendAnalysis(daily2annual(pan_obs_gapfill[istation, :])).KendallTau() for istation in xrange(0, len(r2_60))])
	# pan_mod = load('%s/pan_mod_good_stations' %workspace)
	# panmod = []
	# for ibasin in xrange(0, 10):
	# 	for istation in r2_60_stations[ibasin]:
	# 		panmod.append(pan_mod[istation, :].flatten() * stats_basin[ibasin, 0]) # + stats_basin[ibasin, 1])
	# panmod = vstack(panmod)
	# panmod_p = vstack([TrendAnalysis(daily2annual(panmod[istation, :])).KendallTau() for istation in xrange(0, len(r2_60))])
	# panmod_slp = vstack([TrendAnalysis(daily2annual(panmod[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, len(r2_60))])

	# significance
	# mk_sign = np.where((panobs_p<=0.05)&(panmod_p<=0.05))[0]
	# mk_sign = np.where(panobs_p<=0.05)[0]
	# trend consistency
	# mk_same_sign = np.where(((panobs_slp<0)&(panmod_slp<0))|((panobs_slp>0)&(panmod_slp>0)))[0]
	# mk = intersect1d(mk_significant, mk_same_sign)
	# print r2_60[mk]
	# ##----------------------------------Figure Plotting-----------------------------------------

	# Plotting.Mapshow(panobs_slp[panobs_p>0.05], lons[panobs_p>0.05], lats[panobs_p>0.05], 50, 0.7, min, max, plt.cm.seismic, 'Observed pan evaporation %s' %title, 'p>0.05', unit, figdir, 'pan_obs_%s_r2_60.png' % figname)
	# plt.hold(True)
	# Plotting.Mapshow(panobs_slp[panobs_p<=0.05], lons[panobs_p<=0.05], lats[panobs_p<=0.05], 150, 0.7, min, max, plt.cm.seismic, 'Observed pan evaporation %s' %title, 'p<=0.05', unit, figdir, 'pan_obs_%s_r2_60.png' % figname)

	# Plotting.Mapshow(panmod_slp[panmod_p>0.05], lons[panmod_p>0.05], lats[panmod_p>0.05], 50, 0.7, min, max, plt.cm.seismic, 'Modelled pan evaporation %s' %title, 'p>0.05', unit, figdir, 'pan_obs_%s_r2_60.png' % figname)
	# plt.hold(True)
	# Plotting.Mapshow(panmod_slp[panmod_p<=0.05], lons[panmod_p<=0.05], lats[panmod_p<=0.05], 150, 0.7, min, max, plt.cm.seismic, 'Modelled pan evaporation %s' %title, 'p<=0.05', unit, figdir, 'pan_obs_%s_r2_60.png' % figname)
	#
	# plt.legend(loc=3, prop={'size': 18})
	# plt.show()

	## 3.5 trend significance
	# min = 0
	# max = 1
	# title = 'trend significance'
	# unit = '(mm/yr/yr)'
	# figname = 'pvalue'
	# # pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	# # pan_obs_gapfill = vstack([TrendAnalysis(daily2annual(pan_obs_gapfill[istation, :])).KendallTau() for istation in xrange(0, len(lons))])
	#
	# pan_mod = load('%s/pan_mod_good_stations' %workspace)
	# pan_mod = vstack([TrendAnalysis(daily2annual(pan_mod[istation, :])).KendallTau() for istation in xrange(0, len(lons))])
	# #
	# # ##----------------------------------Figure Plotting-----------------------------------------
	# Plotting.Mapshow('white', lons, lats, 50, min, max, plt.cm.seismic, 'Modelled pan evaporation %s' %title, None, unit, figdir, 'pan_mod_%s.png' % figname)
	# plt.hold(True)
	# Plotting.Mapshow('black', lons[pan_mod<0.1], lats[pan_mod<0.1], 50, min, max, plt.cm.seismic, 'Modelled pan evaporation %s' %title, None, unit, figdir, 'pan_mod_%s.png' % figname)
	# savefig('%s/pan_mod_%s.png' %(figdir, figname))


	# 4. trend attribution
	# min = -10
	# max = 10
	unit = ''
	# figname = 'trend_att_'
	# vars = ['Tair', 'Rn', 'Wind', 'Ea']
	att_stations = load('%s/att_stations' %workspace)
	#
	# # ##----------------------------------Figure Plotting-----------------------------------------
	# [Plotting.Mapshow(att_stations[:, i]*365, lons, lats, 50, 1, min, max, plt.cm.seismic, 'Epan trend attribution to %s' %v, None, unit, figdir, 'pan_mod_trend_att_%s.png' % v) for i, v in enumerate(vars)]

	# # T E+U R
	att_rgb = np.empty((att_stations.shape[0], 3))
	att_rgb[:, 0] = att_stations[:, 0] + att_stations[:, 3]
	att_rgb[:, 1] = att_stations[:, 2]
	att_rgb[:, 2] = att_stations[:, 1]

	attmag = abs(att_rgb)
	maxatt = np.max(abs(att_rgb), axis=1).reshape(att_stations.shape[0], 1)
	att_rgb_norm = attmag / maxatt

	# condition
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)[r2_60,:]
	panobs_slp = vstack([TrendAnalysis(daily2annual(pan_obs_gapfill[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, len(r2_60))])

	pan_mod = load('%s/pan_mod_good_stations' %workspace)
	panmod = []
	for ibasin in xrange(0, 10):
		for istation in r2_60_stations[ibasin]:
			panmod.append(pan_mod[istation, :].flatten() * stats_basin[ibasin, 0]) # + stats_basin[ibasin, 1])
	panmod = vstack(panmod)
	panmod_p = vstack([TrendAnalysis(daily2annual(panmod[istation, :])).KendallTau() for istation in xrange(0, len(r2_60))])
	panmod_slp = vstack([TrendAnalysis(daily2annual(panmod[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, len(r2_60))])
	# trend consistency
	mk_same_sign = np.where(((panobs_slp<0)&(panmod_slp<0))|((panobs_slp>0)&(panmod_slp>0)))[0]
	index_insign = np.where(panmod_p>0.05)[0]
	index_sign = np.where(panmod_p<=0.05)[0]
	mk = intersect1d(mk_same_sign,index_sign)

	# mk_up = intersect1d(np.where((panobs_slp>0)&(panmod_slp>0))[0], index_sign)
	# mk_dw = intersect1d(np.where((panobs_slp<0)&(panmod_slp<0))[0], index_sign)

	mk_up = intersect1d(np.where((panobs_slp>0)&(panmod_slp<0))[0], index_sign)
	mk_dw = intersect1d(np.where((panobs_slp<0)&(panmod_slp>0))[0], index_sign)


	# ##----------------------------------Figure Plotting-----------------------------------------
	## p value as big and small
	# Plotting.Mapshow_RGB(att_rgb_norm[index_insign, :], lons[index_insign], lats[index_insign], 50, 0.4, 'Epan trend attribution : 1961-2001', 'p>0.05', unit)
	# plt.hold(True)
	# Plotting.Mapshow_RGB(att_rgb_norm[index_sign, :], lons[index_sign], lats[index_sign], 150, 0.6, 'Epan trend attribution : 1961-2001', 'p<=0.05', unit)
	## only significant point
	# Plotting.Mapshow_RGB(att_rgb_norm[mk, :], lons[mk], lats[mk], 150, 0.6, 'Epan trend attribution : 1961-2001', 'p<=0.05', unit)
	## significant with up and down triangle
	Plotting.Mapshow_RGB(att_rgb_norm[mk_up, :], lons[mk_up], lats[mk_up], 180, 0.8, '^', 'Epan trend attribution : 1961-2001', 'increase', unit)
	plt.hold(True)
	Plotting.Mapshow_RGB(att_rgb_norm[mk_dw, :], lons[mk_dw], lats[mk_dw], 180, 0.8, 'v', 'Epan trend attribution : 1961-2001', 'decrease', unit)
	plt.legend(loc=4, prop={'size': 14})
	plt.show()
	# savefig('%s/pan_mod_trend_att_%s.png' %(figdir, figname))
	return

Plot_Station_Data()
exit()
def Barplot_trend_attribution_basin():

	## Read trend attribution data
	vars = ['Tair', 'Rn', 'Wind', 'Ea']
	stats_basin = load('%s/calibration_regressionfit_basin_u56' %workspace)
	r2_60_stations = pickle.load(open('r2_60_stations_u56','rb'))
	r2_60 = hstack(r2_60_stations)
	relatetable = load('%s/relatetable' %workspace)
	att_stations = load('%s/att_stations' %workspace)

	# condition
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)[r2_60,:]
	panobs_slp = vstack([TrendAnalysis(daily2annual(pan_obs_gapfill[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, len(r2_60))])

	pan_mod = load('%s/pan_mod_good_stations' %workspace)
	panmod = []
	for ibasin in xrange(0, 10):
		for istation in r2_60_stations[ibasin]:
			panmod.append(pan_mod[istation, :].flatten() * stats_basin[ibasin, 0]) # + stats_basin[ibasin, 1])
	panmod = vstack(panmod)
	panmod_p = vstack([TrendAnalysis(daily2annual(panmod[istation, :])).KendallTau() for istation in xrange(0, len(r2_60))])
	panmod_slp = vstack([TrendAnalysis(daily2annual(panmod[istation, :])).ThSenTrend() * 365.0 for istation in xrange(0, len(r2_60))])
	# trend consistency
	mk_same_sign = np.where(((panobs_slp<0)&(panmod_slp<0))|((panobs_slp>0)&(panmod_slp>0)))[0]
	index_insign = np.where(panmod_p>0.05)[0]
	index_sign = np.where(panmod_p<=0.05)[0]
	mk = intersect1d(mk_same_sign,index_sign)
	mask = empty((att_stations.shape[0], 4))
	mask.fill(np.nan)
	mask[mk,:] = 1
	att_stations = att_stations*mask

	ist = 0
	att_plot = []
	for ibasin in xrange(0, 10):
		att_basin = []
		for istation in relatetable[r2_60_stations[ibasin], 1]:
			att_basin.append(att_stations[ist, :])
			ist = ist + 1
		att_plot.append(nanmean(att_basin, axis=0))
	att_plot = vstack(att_plot) * 365

	## Initialize
	col = ['red', 'blue', 'DarkTurquoise', 'green']
	stackbase_positive = array([0.0]*10)
	stackbase_negative = array([0.0]*10)

	## DRAW FIGURE---------------
	fig, ax = plt.subplots(1, 1, figsize=(16, 7))
	ax.yaxis.grid(True)
	width = 0.6
	bar_positions = np.arange(10)
	for i in xrange(0, 4):
		stackbase_positive = Plotting.Stacked_Bar(ax, att_plot[:, i], width, bar_positions, stackbase_positive, col[i], vars[i], 'positive',0.7)
		stackbase_negative = Plotting.Stacked_Bar(ax, att_plot[:, i], width, bar_positions, stackbase_negative, col[i], vars[i], 'negative',0.7)


	# Set up xaxis
	plt.xlim([0, 10.2])
	ax.set_xticks(bar_positions + width)  # add xtick
	ax.set_xticklabels(basinlongs, fontsize=24, fontweight='bold', rotation=30)
	ax.tick_params(axis='x', length=0) # leave space between pad pad = 10
	# Set up yaxis
	ax.set_ylabel('[mm/yr/yr]', fontsize=24, fontweight='bold')
	ax.tick_params(axis='y', labelsize=24)

	# make the top and right axes invisible
	ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.spines['bottom'].set_visible(False)
	ax.get_xaxis().tick_bottom(); ax.get_yaxis().tick_left()

	# add the direction text
	# ymin, ymax = ax.get_ylim()
	# ax.text(0.3, ymax-0.9*(ymax-ymin), "North", fontsize=25, fontweight='bold')
	# ax.text(4.4, ymax-0.9*(ymax-ymin), "South", fontsize=25, fontweight='bold')

	# add legend below
	lgd = plt.legend(loc='upper right', fontsize=20)
	savefig('%s/mod_trend_attribution_basin_barplot.png' %figdir, bbox_extra_artists=(lgd, ), bbox_inches='tight')
	plt.show()

	return
# Barplot_trend_attribution_basin()





