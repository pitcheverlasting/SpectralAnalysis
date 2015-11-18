__author__ = 'lpeng'

from pylab import *
from scipy import signal
import scipy
import math

def Power_Spectrum(ts, fs, dt):

	"""Compute power spectral density (using psd welch)"""
	"ts: time series of measurement values"
	"fs: sampling frequency"
	"dt: detrend = None/'linear'/'constant'"
	"scal = 'density', normalized with units of V**2/hz"
	"scal = 'spectrum', units of V**2"

	ts[np.isnan(ts)==True] = nanmean(ts)		# Gapfill the missing value
	nwindow = 9          						# Number of windows to smooth data
	length = math.floor(len(ts)/nwindow)  		# Length calculated by deviding the window
	nwindow_fl = math.floor(log2(length))		# Number of windows with floor window length
	window = int(2 ** nwindow_fl)   			# segment_length

	############################### nfft: Number of DFT points ###############################
	# The default nfft is the greater of 256 or the next power of 2 greater than the length of the segments.
	# The relation between segment length and the Number of DFT points: segment_length <= nfft
	# If nfft is greater than the segment length, the data is zero-padded. If nfft is less than the segment length, the segment is wrapped using datawrap to make the length equal to nfft
	##########################################################################################

	###############################   welch's method for PSD   ###############################
	# MATLAB function: pwelch(ts, segment_length, Number of overlapped window noverlap, Number of DFT points, fs,'onesided')
	# scipy method: pwelch(ts, fs, nperseg, noverlap, nfft, detrend, scaling)
	##########################################################################################

	f, Pxx_den = signal.welch(ts, fs, nperseg=window, detrend=dt)

	return [f, Pxx_den]


def CrossPowerSpectrum(ts1, ts2, fs, dt):

	"""Compute Cross spectral density (using csd)"""

	# Gapfill the missing value
	ts1[np.isnan(ts1)==True | np.isnan(ts2)==True] = nanmean(ts1)
	ts2[np.isnan(ts1)==True or np.isnan(ts2)==True] = nanmean(ts2)
	nwindow = 9          						# Number of windows to smooth data
	length = math.floor(len(ts1)/nwindow)  		# Length calculated by deviding the window
	nwindow_fl = math.floor(log2(length))		# Number of windows with floor window length
	window = int(2 ** nwindow_fl)   			# segment_length

	f, Pxy_den = signal.csd(ts1, ts2, fs, nperseg=window, detrend=dt)

	return [f, Pxy_den]


def Coherence(ts1, ts2, fs, dt):

	"""Compute magnitude squared coherence (using coherence)"""

	# Gapfill the missing value
	mk1 = np.isnan(ts1)
	mk2 = np.isnan(ts2)
	ts1[mk1==True] = nanmean(ts1)
	ts1[mk2==True] = nanmean(ts1)
	ts2[mk1==True] = nanmean(ts2)
	ts2[mk2==True] = nanmean(ts2)

	nwindow = 9          						# Number of windows to smooth data
	length = math.floor(len(ts1)/nwindow)  		# Length calculated by deviding the window
	nwindow_fl = math.floor(log2(length))		# Number of windows with floor window length
	window = int(2 ** nwindow_fl)   			# segment_length

	f, Cxy = signal.coherence(ts1, ts2, fs, nperseg=window, detrend=dt)

	return [f, Cxy]

def WaveletCoherence():

	return