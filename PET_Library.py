__author__ = 'lpeng'
"""Forcing Preprocess library"""

import pandas as pd
from pylab import *

class Data:

	# Initialization
	def __init__(self, grid_data):

		# Read data from the binary file
		# forcing = ('tair', 'tdew', 'pmsl', 'wdir', 'wspd', 'skyc', 'prc1', 'prc6', 'pmst', 'shum')
		# forname = ('Air Temperature (Degrees Celsius, x10)', 'Wind Speed Rate (meters per second, x10)', 'Station Level Pressure (Hectopascals, x10)', 'Specific Humidity - (g/kg, x100)')
		INPUT = self.Read_Grid_Data(grid_data)

		# Define the incoming grid data variables
		self.Tmax = INPUT['Tmax']
		self.Tmin = INPUT['Tmin']
		self.RH = INPUT['RH']
		# self.Wind = self.Convert_Unit_Wind(INPUT['wnd10m'])
		self.lat = INPUT['lat']
		self.doy = INPUT['doy']
		self.gslen = INPUT['gslen']
		# self.Tair = self.Convert_Unit_Temp(INPUT['Tair'])
		# self.Pres = self.Convert_Unit_Pres(INPUT['Pres'])
		# PenPan package: P = 101.3 * ((293 - 0.0065 * Elev)/293) ** 5.26
		self.SHum = INPUT['SHum']  # (kg kg^-1)
		self.Wind = INPUT['Wind']  # Wind: m/s
		self.CF = INPUT['tc']
		self.Albedo = INPUT['Albedo']

		# Calculate some variables
		self.Calculate_Rnet()
		self.Calculate_Saturated_Vapor_Pressure()
		self.Calculate_Humidity()  # self.e
		self.Calculate_Slope_Saturation_Vapor_Pressure_Curve()
		self.Calculate_Gamma()
		self.Calculate_VPD()
		self.Penman()


	def Read_Grid_Data(self, grid_data):

		varname = ('LWin', 'SWin', 'Tair', 'Wind', 'Pres', 'SHum', 'Albedo')
		INPUT = {varname[i]: grid_data[i, :, :] for i in range(0, len(varname))}

		return INPUT

	def Calculate_Saturated_Vapor_Pressure(self):

		self.estar = 0.6108 * np.exp((17.27 * self.Tair) / (237.3 + self.Tair))

		return

	def Calculate_Tmean(self):

		self.Tmean = (self.Tmax + self.Tmin) / 2.0

		return

	def Calculate_Mean_Saturated_Vapor_Pressure(self):

		estar_max = 0.6108 * np.exp((17.27 * self.Tmax) / (237.3 + self.Tmax))
		estar_min = 0.6108 * np.exp((17.27 * self.Tmin) / (237.3 + self.Tmin))
		self.estar = (estar_max + estar_min) / 2

		return

	def Calculate_Humidity(self):

		mixing = self.SHum / (1 - self.SHum)
		self.e = mixing * self.Pres / (0.622 + mixing)
		self.e = self.estar * (self.estar<self.e) + self.e * (self.estar>=self.e)

		return

	def Calculate_Slope_Saturation_Vapor_Pressure_Curve(self):

		# self.DELTA = 4098 * self.estar / (237.3 + self.Tair) ** 2
		# self.DELTA = 4098 * self.estar / (237.3 + self.Tmean) ** 2
		self.DELTA = 4098 * 0.6108 * np.exp((17.27 * self.Tmean) / (237.3 + self.Tmean)) / (237.3 + self.Tmean) ** 2

		return

	def Calculate_Rnet(self):
		"Compute Rnet with longwave downward and shortwave downward radiation data"
		# albedo = 0.23
		# self.Rn = (1 - albedo) * self.SWin + self.LWin - stefan_b * (self.Tair+273.16)**4

		stefan_b = 4.903e-9  # [MJ K-4 m-2 day-1]
		# add mask to long wave radiation
		mask = ones((self.SWin.shape[0], self.SWin.shape[1]))
		mk = np.isnan(self.SWin)
		mask[mk==True] = np.nan
		stefan_b * (self.Tair+273.16)**4 * mask
		self.LWnet = stefan_b * ((self.Tmax+273.16)**4 + (self.Tmin+273.16)**4) / 2.0 * (0.34-0.14 * np.sqrt(self.e)) * (1.35 * self.CF - 0.35)
		self.Rn = (1 - self.Albedo) * self.SWin + self.LWnet

		return

	def Calculate_Rso(self):

		dr = 1 + 0.033 * np.cos(2 * np.pi * self.doy/365.0)  # dr is inverse relative distance Earth-Sun
		d = 0.409 * np.sin((2*np.pi / 365.0) * self.doy - 1.39)   # d is solar declination
		phi =(np.pi/180) * self.lat   # from decimal degree to radians
		omega = np.arccos(-np.tan(phi) * np.tan(d))  # sunset hour angle
		if self.elev:
			z2 = self.elev
		else:
			z2 = 30  # set the station elevation above the sea level as 30m
		Gsc = 0.082   # solar constant = 0.082 [MJ m^-2 min^-1]
		# Ra: extraterrestrial radiation for daily period
		self.Ra = 24 * 60 / np.pi * Gsc*dr*(omega * np.sin(phi) * np.sin(d) + np.cos(phi) * np.cos(d) * np.sin(omega))
		Rso = (0.75 + 2e-5 * z2) * self.Ra
		## this is a temporary function for Rso
		self.Rso = np.tile(Rso, self.nt/self.gslen)

		return

	def Calculate_Rs_pan(self):
		"Compute the total SW irradiance of the pan"
		"Prad : pan evaporation factor, which accounts for the extra direct irradiance intercepted by the walls of the pan when the sun is not directly overhead"
		"f_dir : fraction of fs that is direct"
		"R_s : downward solar irradiance at the surface"
		# if solar == "data":
		# 	R_s = self.SWin
		# elif solar == "sunhours":
		# 	# R_s = (as + bs * (n/N)) * self.Ra
		# else:
		# 	# R_s = (0.85 - 0.047 * Cd) * self.Ra
		# N = 24/pi * w_s


		# ET_type = "Class-A Pan Evaporation"   ("Evaporative surface: ", Surface)
		vars = ['Tmax', 'Tmin', 'RHmax', 'RHmin', 'u2', 'uz', "sunshine hours", "cloud", "monthly precipitation"]

		P_rad = 1.32 + 4 * 10**(-4) * abs(self.lat) + 8 * 10 ** (-5) * self.lat ** 2
		f_dir = -0.11 + 1.31 * self.Rs/self.Ra
		self.Rsp = (f_dir * P_rad + 1.42 * (1 - f_dir) + 0.42 * self.Albedo) * self.Rs

		return

	def Calculate_Rnet_pan(self):

		Ap = 0.14  # Class A pan albedo (Linacre, 1992; Rotstayn, 2006; Roderick, 2007; Yang and Yang, 2011)
		stefan_b = 4.903e-9  # [MJ K-4 m-2 day-1]
		self.LWnet = stefan_b * ((self.Tmax+273.16)**4 + (self.Tmin+273.16)**4) / 2.0 * (0.34-0.14 * np.sqrt(self.e)) * (1.35 * self.CF - 0.35)
		self.SWout = self.Albedo * self.SWin
		self.Rn_pan = (1 - Ap) * self.Rsp - self.LWnet

		return

	def Calculate_Gamma(self):

		cp = 1.013   # Specific heat of moist air at constant pressure [kJ kg-1 C-1]
		self.lv = 2.501 - 0.002361 * self.Tair  # Latent heat vapor
		self.gamma = ((cp * self.Pres) / (0.622 * self.lv)) * 10 ** (-3)
		# PenPan package: gamma = 0.00163 * P/lambda

		return self.gamma

	## Convert unit
	# kelvin to degree
	def Convert_Unit_Temp(self, input):
		data = input - 273.16
		return data

	# W/m2 to MJ/m2/day
	def Convert_Unit_Rad(self, input):
		watt2jule = 10e5/86400.0
		data = input / float(watt2jule)
		return data

	# pressure: pa to kpa

	def Convert_Unit_Pres(self, input):
		data = input / 1000.0
		return data

	# 10m wind to 2m wind
	def Convert_Unit_Wind(self, input):

		zh = 10  # 10 m wind field
		data = (4.87 / (np.log(67.8 * zh - 5.42))) * input

		return data

	## Calculation for each components

	def Calculate_VPD(self):

		self.vpd = self.estar - self.e

		return self.vpd

	# Reference-surface Models
	def Penman(self):

		# Hydrology Book
		# all data are half-hourly time step
		PET_R = (self.DELTA / (self.DELTA + self.gamma)) * self.Rn / self.lv
		PET_A = (self.gamma / (self.DELTA + self.gamma)) * ((6.43 * (1 + 0.536 * self.Wind) * self.vpd) / self.lv)
		self.PET_Penman = PET_R + PET_A

		return

	def Penpan(self, solar):


		"alpha: Please use a numeric value for the alpha (albedo of evaporative surface)"
		"Solar radiation data have been used directly for calculating evapotranspiration"
		"Sunshine hour data have been used for calculating incoming solar radiation"
		"Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
		"Monthly precipitation data have been used for calculating incoming solar radiation"


		# vabar = (vs_Tmin * RHmax/100 + vs_Tmax * RHmin/100)/2
		#
		coeff_hq = 2.4
		f_pan_u = 1.39e-8 * (1 + 1.35 * self.Wind)
		PET_R = self.DELTA/(self.DELTA + coeff_hq * self.gamma) * self.Rn_pan / self.lv
		PET_A = coeff_hq * self.gamma/(self.DELTA + coeff_hq * self.gamma) * f_pan_u * self.vpd
		self.penpan = PET_R + PET_A
		# if overest == TRUE:
		# 	Epenpan.Daily = Epenpan.Daily/1.078

    	return