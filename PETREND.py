__author__ = 'lpeng'
"""PET trend analysis library"""

from pylab import *
import numpy as np
import pandas as pd
from scipy import stats

class TrendAnalysis(object):

	def __init__(self, input):

		self.tstep = len(input)
		self.INPUT = input

	def ThSenTrend(self):

		""" Use Theil-Sen slope to calculate the median (nonparametric) trend """

		trend = stats.mstats.theilslopes(self.INPUT, alpha=0.95)[0]  # Recall that mstats use alpha=0.95 but stats use alpha=0.05

		return trend

	def KendallTau(self):

		""" Use Kendall robust tau test to calculate the trend significance """

		corr, p = stats.mstats.kendalltau(np.arange(1, self.tstep+1), self.INPUT)

		return p


class TrendAttribution(object):

	def __init__(self, forave, fortrend):

		self.AVE = forave
		self.TREND = fortrend
		
	def Calculate_Attribution(self):

		# Tair_ave = self.AVE['Tair_ave']
		gamma_ave = self.AVE['gamma_ave']
		delta_ave = self.AVE['delta_ave']
		vpd_ave = self.AVE['vpd_ave']
		wind_ave = self.AVE['wind_ave']
		lhv_ave = self.AVE['lhv_ave']
		Rn_ave = self.AVE['Rn_ave']
		PDRnl_ave = self.AVE['PDRnl_ave']

		trend_Rnet = self.TREND['trend_Rn']
		trend_Delta = self.TREND['trend_Delta']
		trend_Wind = self.TREND['trend_Wind']
		trend_Esat = self.TREND['trend_Esat']
		trend_Eact = self.TREND['trend_Eact']
		trend_Tair = self.TREND['trend_Tair']

		# calculate all partial derivatives for each components
		coeff_hq = 5
		# f_pan_u = 2.6 * (1 + 0.536 * Wind)
		# Df_pan_u = 2.6 * 0.536
		# 0.35*(1+9.8e-3 * self.Wind)
		# 1.39e-8 * 1.35
		# 1.313 + 1.381 * self.Wind
		Df_pan_u = 1.381
		# Epan = Epan,Rn + Epan,Aero
		# A. derivative for Rnet, in Epan,Rn
		PD_Rn = (delta_ave / (delta_ave + coeff_hq * gamma_ave)) * trend_Rnet / lhv_ave
		# Net Radiation Rn:  (delta_ave, gamma_ave, trend_Rnet)
		dRn = PD_Rn

		# B. derivative for U2 (wind), in Epan,Aero
		PD_u = ((Df_pan_u * coeff_hq * gamma_ave * vpd_ave) / (delta_ave + coeff_hq * gamma_ave)) * trend_Wind
		# Wind speed U:  (delta_ave, gamma_ave, vpd_ave, lhv_ave, trend_Wind)
		dU = PD_u

		# C. derivative for e_a (humidity), in Epan,Aero
		PD_ea = ((- coeff_hq * gamma_ave * (1.313 + 1.381 * wind_ave)) / (delta_ave + coeff_hq * gamma_ave)) * trend_Eact
		# Humidity Ea: (gamma_ave, wind_ave, delta_ave, lhv_ave, trend_Eact)
		dEa = PD_ea

		# D. derivative for temperature
		# D1. derivative for temperature, in Epan, Rn: delta/(delta+gamma) part
		PDR_Ddelta = ((Rn_ave * coeff_hq * gamma_ave) / (delta_ave + coeff_hq * gamma_ave) ** 2) * trend_Delta / lhv_ave
		dRn_dTair = PDRnl_ave
		PDR_DRn = (delta_ave / (delta_ave + coeff_hq * gamma_ave)) * dRn_dTair * trend_Tair / lhv_ave
		# D2. derivative for temperature, in Epan, Aero: gamma/(delta+gamma) part
		PDA_Ddelta = (- coeff_hq * gamma_ave * vpd_ave * (1.313 + 1.381 * wind_ave)) / ((delta_ave + coeff_hq * gamma_ave) ** 2) * trend_Delta
		PDA_Des = (coeff_hq * gamma_ave * (1.313 + 1.381 * wind_ave)) / ((delta_ave + coeff_hq * gamma_ave)) * trend_Esat
		# Temperature Ta: (Rn_ave, gamma_ave, delta_ave, trend_Delta, trend_Tair, vpd_ave, wind_ave, lhv_ave, trend_Esat)
		dTa = PDR_Ddelta + PDR_DRn + PDA_Ddelta + PDA_Des


		return dTa, dRn, dU, dEa
