import astropy.constants as ac
import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as au
import astropy.constants as ac
from astropy.io import fits
import math 
import scipy
from scipy.interpolate import interp2d
from matplotlib.pyplot import cm
import matplotlib
import scipy
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import struct
import lmfit
from lmfit import Model, Parameter, report_fit

class binary_file:
	def __init__(self, filename):
		self.filename = filename

	def exp_list(self, start, stop, total_number, base):
		"""Calculates exponential distribution between two values.

		Args:
			start (float): starting value
			stop (float): last value
			total_number (int): total amount of distributed values.
			base (float): distribution factor


		Returns:
			number_list (list): Distributed numbers
		"""

		number_list = np.zeros(total_number + 1)
		number_list[0] = start
		number_list[total_number] = stop

		if base > 1:
			dx = (stop - start) * (base - 1.0) / (pow(base, total_number) - 1)

			for i_x in range(0, total_number):
				number_list[i_x] = start + dx * \
					(pow(base, i_x) - 1) / (base - 1.0)
		else:
			raise ValueError('only positive exp bases are allowed!')
		return number_list


	def radial_list(self, dictionary):
		"""
		Args:
			Dictionary of the grid

		Returns:
			Position of the center of the cell in the radial direction
		"""
		radial_list = self.exp_list(dictionary['Rmin'], dictionary['Rmax'], dictionary['Nr'], dictionary['fr'])/au.au.to('m')
		length_cell = radial_list[1:] - radial_list[0:len(radial_list)-1]
		radial_list = radial_list[0:len(radial_list)-1] + length_cell/2
		return radial_list
	
	def parameter_type(self, ID, idx, parameter):
		if ID == 0:
			parameter['n_H'] = idx
		elif ((ID == 2) and (idx == 1)):
			parameter['T_d'] = idx
		elif ID == 4:
			parameter['B_x'] = idx
		elif ID == 5:
			parameter['B_y'] = idx
		elif ID == 6:
			parameter['B_z'] = idx
		
		elif ((ID == 2) and (idx != 1)):
			parameter['T_d_a'].append(idx)
		elif ID == 47:
			parameter['abs_ini'].append(idx)
		elif ID == 13:
			parameter['a_align'] = idx
		elif ID == 28:
			parameter['m_H'] = idx
		elif ID == 30:
			parameter['ux'].append(idx)
		elif ID == 31:
			parameter['uy'].append(idx)
		elif ID == 32:
			parameter['uz'].append(idx)
		elif ID == 33:
			parameter['urad'].append(idx)
		elif ID == 34:
			parameter['cos(psi)'] = idx		
		elif ID == 35:
			parameter['gamma'] = idx
		elif ID == 36:
			parameter['adisr'] = idx
		elif ID == 37:
			parameter['eta'] = idx
		elif ID == 38:
			parameter['adisr_max'] = idx
		elif ID == 39:
			parameter['a_min_aJ_lowJ'] = idx
		elif ID == 40:
			parameter['a_max_aJ_lowJ'] = idx
		elif ID == 41:
			parameter['a_min_aJ_highJ'] = idx
		elif ID == 42:
			parameter['a_max_aJ_highJ'] = idx
		elif ID == 43:
			parameter['a_min_JB_DG_50'] = idx
		elif ID == 44:
			parameter['a_max_JB_DG_50'] = idx		
		elif ID == 45:
			parameter['a_min_JB_DG_100'] = idx
		elif ID == 46:
			parameter['a_max_JB_DG_100'] = idx
		elif ID == 48:
			parameter['a_max_JB_Lar'] = idx	
		elif ID == 49:
			parameter['akrat_lowJ'] = idx
		elif ID == 50:
			parameter['akrat_highJ'] = idx	
				
	def read_parameter(self, ID_Nphys, parameter):
		parameter['T_d_a'] = []
		parameter['abs_ini'] = []
		parameter['urad'] = []
		parameter['ux'] = []
		parameter['uy'] = []
		parameter['uz'] = []
		for idx, ID in enumerate(ID_Nphys):
			self.parameter_type(ID, idx, parameter)
		return parameter
		
	def read_binary_grid_file(self):
		with open(self.filename, 'rb') as fp:
		 
			# Basic information of grid
			grid_id, = (struct.unpack('<H', fp.read(2)))
			parameter_size, = struct.unpack('<H', fp.read(2))
			ids = struct.unpack('<' + parameter_size*'H', fp.read(2*parameter_size))
			Rmin, = struct.unpack('<d', fp.read(8))
			Rmax, = struct.unpack('<d', fp.read(8))
			Nr, = struct.unpack('<H', fp.read(2))
			Nph, = struct.unpack('<H', fp.read(2))
			Nth, = struct.unpack('<H', fp.read(2))
			fr, = struct.unpack('<d', fp.read(8))
			fph, = struct.unpack('<d', fp.read(8))
			fth, = struct.unpack('<d', fp.read(8))

			# Dictionary of header
			parameter = dict()
			parameter['parameter_size_cell'] = parameter_size
			parameter['Rmin'] = Rmin
			parameter['Rmax'] = Rmax
			parameter['Nr'] = Nr
			parameter['Nph'] = Nph
			parameter['Nth'] = Nth
			parameter['fr'] = fr
			parameter['fph'] = fph
			parameter['fth'] = fth
			
			# Take data in all cells in grid and reshape it into 4D matrix
			data = struct.unpack('<' + parameter_size * Nr * Nth * Nph * 'd', fp.read(8 *  parameter_size * Nr * Nth * Nph))  
			center_point = struct.unpack('<' + parameter_size * 'd', fp.read(8 * parameter_size))
			matrix = np.reshape(data, (Nr, Nph, Nth, parameter_size))
			
			# Dictionary for parameter inside file
			parameter = self.read_parameter(ids, parameter)
			if len(parameter['urad']) != 0:
				parameter['U'] = self.radiation_field(parameter, matrix)
		return parameter, matrix	

	def radiation_field(self, parameter, matrix, nr_wave = 100):
		wavelength = self.wavelength_list() 
		Wave = np.zeros(nr_wave)
		Ulambda = np.zeros((parameter['Nr'], parameter['Nph'], parameter['Nth'], nr_wave))
		for i in range(nr_wave):
			Wave[i] = wavelength[i]
			Ulambda[:, :, :, i] = matrix[:, :, :, parameter['urad'][i]]/ac.c.value
		U = np.trapz(Ulambda, Wave, axis = 3)
		return U/8.64e-14		 
		 
		 
		 
	def wavelength_list(self):
		"""
		Returns: list of wavelengths
		"""
		wave_list = np.genfromtxt(os.path.join('/Users/chaugiang/Dropbox/POLARIS-/input/wave.dat'))
		return wave_list
		
	def radius_list(self):
		"""
		Returns: list of wavelengths
		"""
		radius_list = np.genfromtxt(os.path.join('/Users/chaugiang/Dropbox/POLARIS-/input/radius.dat'))
		return radius_list
