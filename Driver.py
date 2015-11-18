__author__ = 'lpeng'

from pylab import *
import scipy.io
import pandas as pd
import IO, FFT, Plotting

# #=============================path================================
datadir = '/home/water5/lpeng/script/PROJECT/pan_evaporation_china/original'
workspace = '/home/water5/lpeng/script/PROJECT/pan_evaporation_china/workspace'
##=============================variable============================
vars = ['time', 'p', 'tavg', 'tmin', 'tmax', 'ea', 'rh', 'tc', 'lc', 'wind', 'ts', 'sun', 'rain', 'pan', 'vpd',  'estavg', 'rh_test']
varlongs = ['time', 'daily averaged pressure [hpa]', 'daily averaged air temperature [degree]', 'daily minimum air temperature [degree]', 'daily maximum air temperature [degree]', 'daily averaged vapor pressure [hpa]', 'daily averaged relative humidity [%]', 'daily averaged total cloud cover', 'daily averaged low cloud cover', 'daily averaged wind speed [m/s]', 'daily averaged surface temperature(at 0cm level) [degree]', 'daily averaged sunshine hour [hr]', 'daily averaged rainfall [mm]', 'daily averaged pan evaporation [mm]', 'daily averaged vapor pressure deficit [hpa]', 'daily saturated vapor pressure using daily averaged air temperature [hpa]', 'daily averaged relative humidity calculated [%]']
# Note that original data is multiplied by 10 from the normal unit. This process data devide 10 so that it becomes the original unit
# time 1,2,3 year, month, day
# for coherece analysis
variables = ['tavg', 'sun', 'wind', 'ea', 'lc', 'ts', 'vpd', 'rain']

basinlongs=['Songhuajiang','Liaohe','Northwestern','Haihe','Yellow', 'Yangtze','Huaihe','Southeastern','Southwestern','Pearl']
# geoinfo = load('%s/station_geoinfo' %workspace)
station_number = load('%s/basin_station_number' %workspace)

# Dimensions
styr = 1961
edyr = 2006
stdy = datetime.datetime(styr, 1, 1)
eddy = datetime.datetime(edyr, 12, 31)
dates = pd.date_range(stdy, eddy, freq='D')
tstep = len(dates)
nf = 513
nbasin = 10
nvar = 8
sampling_frequency = 1/(24.0 * 3600.0)  # unit: per day


def Coherence_Station():

	for ibasin in xrange(1, 2): #11):
		for istation in xrange(0, 1):  #station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin))
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [data[v][0, 0].flatten() for v in variables]
			pan = data['pan'][0, 0].flatten()
			# result = FFT.Power_Spectrum(pan, sampling_frequency, 'linear')
			coh = [FFT.Coherence(var, pan, sampling_frequency, 'linear')[1] for var in input]
			freq = FFT.Coherence(input[0], pan, sampling_frequency, 'linear')[0]
			Plotting.CoherenceStationPlot(coh, sampling_frequency, freq)

	return

## method1: average the coherence spectrum
def Coherence_Station2Basin():

	fig = plt.figure(figsize=(12, 30))
	cohere = []

	for ibasin in xrange(0, 10):
		cohere_basin = zeros((station_number[ibasin], nvar, nf))

		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			## Quality check
			input = [np.isnan(nanmean(data[v][0, istation])) for v in variables]
			flag = []
			flag = [1 for i in input if i == True]
			if len(flag) > 0:
				print "Ignore %s station %s!" % (basinlongs[ibasin], istation)
				continue

			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [data[v][0, istation].flatten() for v in variables]
			pan = data['pan'][0, istation].flatten()
			# Compute the coherence
			freq = FFT.Coherence(input[0], pan, sampling_frequency, 'linear')[0]
			cohere_basin[istation, :, :] = vstack([FFT.Coherence(v, pan, sampling_frequency, 'linear')[1] for v in input])

		# store basin average
		cohere.append(mean(cohere_basin, axis=0).reshape(1, nvar, nf))

		## DRAW FIGURE---------------
	 	ax = fig.add_subplot(4, 3, ibasin+1)
		Plotting.CoherenceBasinPlot(ax, mean(cohere_basin, axis=0), sampling_frequency, freq, basinlongs[ibasin])

	# for national average
	ax = fig.add_subplot(4, 3, 11)
	Plotting.CoherenceBasinPlot(ax, mean(vstack(cohere), axis=0), sampling_frequency, freq, 'Average')
	ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left')
	# fig.tight_layout()
	plt.show()

	return

## method2: average the time series first and then do the spectrum
def Coherence_average2Basin():

	fig = plt.figure(figsize=(12, 30))
	input = []; pan = []; cohere = []

	for ibasin in xrange(0, 10):
		input_basin = zeros((station_number[ibasin], nvar, tstep))
		pan_basin = zeros((station_number[ibasin], tstep))
		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			## Quality check
			input_mask = [np.isnan(nanmean(data[v][0, istation])) for v in variables]
			flag = []
			flag = [1 for i in input_mask if i == True] # which means one variable is missing
			if len(flag) > 0:
				print "Ignore %s station %s!" % (basinlongs[ibasin], istation)
				continue

			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input_basin[istation, :, :] = vstack([data[v][0, istation].flatten() for v in variables])
			pan_basin[istation, :] = data['pan'][0, istation].flatten()

		input.append(mean(input_basin, axis=0).reshape(1, nvar, tstep))
		pan.append(mean(pan_basin, axis=0))
		del input_basin, pan_basin
		# Compute the coherence
		freq = FFT.Coherence(input[ibasin][0, 0, :], pan[ibasin], sampling_frequency, 'linear')[0]
		coh = vstack([FFT.Coherence(input[ibasin][0, v, :], pan[ibasin], sampling_frequency, 'linear')[1] for v in xrange(0, nvar)])  # nvar, nf
		# store basin average
		cohere.append(coh.reshape(1, nvar, nf))

		## DRAW FIGURE---------------
	 	ax = fig.add_subplot(4, 3, ibasin+1)
		Plotting.CoherenceBasinPlot(ax, coh, sampling_frequency, freq, basinlongs[ibasin])

	# for national average
	# method 1: calculate national average first then coherence
	input_all = mean(vstack(input), axis=0) # nvar, tstep
	pan_all = mean(vstack(pan), axis=0).flatten() # 1, tstep
	coh1 = vstack([FFT.Coherence(input_all[v, :].flatten(), pan_all, sampling_frequency, 'linear')[1] for v in xrange(0, nvar)])
	ax = fig.add_subplot(4, 3, 11)
	Plotting.CoherenceBasinPlot(ax, coh1, sampling_frequency, freq, 'Average1')
	# method 2: average the basin coherence to get the national one
	coh2 = mean(vstack(cohere), axis=0)
	ax = fig.add_subplot(4, 3, 12)
	Plotting.CoherenceBasinPlot(ax, coh2, sampling_frequency, freq, 'Average2')
	# ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left')
	# fig.tight_layout()
	plt.show()

	return

##====================================================================
# Run the step functions
##====================================================================
# Coherence_Station()

# Coherence_Station2Basin()
Coherence_average2Basin()