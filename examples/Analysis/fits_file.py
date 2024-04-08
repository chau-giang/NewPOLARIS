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

class plot_2d:

	def __init__(self, filename):
		"""
		type_simulation: [str] "mcmc" or "ray_tracing"
		type_fits_file: [str] "input" or "output"
		"""

		#if len(inspect.getargspec(plot_3d)[0]) != 3:
		#	print('This function includes 2 arguments: filename and type_simulation')
		#	print('For type_simulation: please enter "mcmc" or "ray_tracing"')
		self.filename = filename 
		

	def take_data(self, ID, header, data, dictionary):
		# Gas density, magnetic fields
		#-------------------------------------------------------------------------------
		if 'gas_mass_density [kg/m^3]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:		
			dictionary['mH'] = data[ID, :, :, :] #um
			avg_atomic_mass = 2
			dictionary['nH'] = dictionary['mH']*1e-6/(avg_atomic_mass * ac.m_p.value)
		elif 'mag_total [T]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['B'] = data[ID, :, :, :] #um		 	 		
		elif 'mag_x [T]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['B_x'] = data[ID, :, :, :]*1e10 #um		 
		elif 'mag_y [T]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['B_y'] = data[ID, :, :, :]*1e10 #um		 
		elif 'mag_z [T]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['B_z'] = data[ID, :, :, :]*1e10 #um		 
	 						 
		# Gas temperature, Dust temperature
		#-------------------------------------------------------------------------------
		elif 'gas_temperature [K]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
	 		dictionary['Tg'] = data[ID, :, :, :] #um
		elif 'dust_temperature [K]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['Td'] = data[ID, :, :, :] #um
	
	
		# Radiation field information
		#-------------------------------------------------------------------------------
		elif 'u_rad/u_isrf' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['U'] = data[ID, :, :, :]
		elif 'avg. RAT cos(theta)'in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['theta'] = np.arccos(data[ID, :, :, :])*180/math.pi
		elif 'avg. RAT aniso. (gamma)'in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['gamma'] = data[ID, :, :, :]
	
	
		# Grain alignment by RATs
		#-------------------------------------------------------------------------------
		elif 'rat_aalig [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['align'] = data[ID, :, :, :]*1e6 #um
		elif 'amaxJB_Lar [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['amaxJB_Lar'] = data[ID, :, :, :]*1e6 #um	
	
	
		# Grain alignment by MRAT 
		#-------------------------------------------------------------------------------
		elif 'amin_aJ_lowJ [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['abar_low_min'] = data[ID, :, :, :]*1e6 #um
		elif 'amax_aJ_lowJ [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['abar_low_max'] = data[ID, :, :, :]*1e6 #um
		elif 'amin_aJ_highJ [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['abar_high_min'] = data[ID, :, :, :]*1e6 #um
		elif 'amax_aJ_highJ [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['abar_high_max'] = data[ID, :, :, :]*1e6 #um
		elif 'aminJB_DG_0.5 [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['adg_50_min'] = data[ID, :, :, :]*1e6 #um
		elif 'amaxJB_DG_0.5 [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['adg_50_max'] = data[ID, :, :, :]*1e6 #um
		elif 'aminJB_DG_1[m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['adg_100_min'] = data[ID, :, :, :]*1e6 #um
		elif 'amaxJB_DG_1 [m]' in header['HIERARCH MIDPLANE{:d}'.format(ID+1)]:
			dictionary['adg_100_max'] = data[ID, :, :, :]*1e6 #um
	 	 
		 
	def read_fits(self):
		hdulist = fits.open(self.filename)
		hdu = hdulist[0]

		dat = dict()
		dat['boundary'] = dict()
		dat['boundary']['Nx'] = hdu.header['NAXIS1']
		dat['boundary']['Ny'] = hdu.header['NAXIS2']
		dat['boundary']['Nz'] = hdu.header['NAXIS3']
		dat['boundary']['N_phys'] = hdu.header['NAXIS4']
		
		dat['boundary']['xmin'] = hdu.header['CRVAL1B']
		dat['boundary']['xmax'] = abs(dat['boundary']['xmin'])
		dat['boundary']['dx'] = hdu.header['CDELT1B']

		dat['boundary']['ymin'] = hdu.header['CRVAL2B']
		dat['boundary']['ymax'] = abs(dat['boundary']['ymin'])
		dat['boundary']['dy'] = hdu.header['CDELT2B']
		
		
		dat['data'] = dict()
		for ID in range (dat['boundary']['N_phys']):
			self.take_data(ID, hdu.header, hdu.data, dat['data'])
			
		if "output_midplane.fits" in self.filename:
			# take data from input_midplanes_fits file
			input_fname = self.filename.replace('output','input')
			hdulist = fits.open(input_fname)
			hdu = hdulist[0]
			ID = 0 # ID of gas density
			self.take_data(ID, hdu.header, hdu.data, dat['data'])
			
		return dat

		
def configure_figure(ax):
# Configure figureƒbina	
	for axis in ['top', 'bottom','left','right']:
		ax.spines[axis].set_linewidth(1.5)
		ax.tick_params(bottom=True, top=False, left=True, right=True)
		ax.xaxis.set_tick_params(labelsize=13, width=2)
		ax.xaxis.set_tick_params(which='minor',width=1.15, size=4, direction='in', bottom = True, top = False)
		ax.xaxis.set_tick_params(which='major',width=1.15, size=9, direction ='in', bottom = True, top = False)
		ax.yaxis.set_tick_params(labelsize=13, width=2)
		ax.yaxis.set_tick_params(which='minor',width=1.15, size=4, direction='in', left = True, right = False)
		ax.yaxis.set_tick_params(which='major',width=1.15, size=9, direction='in', left = True, right = False)	 
  
def configure_colorbar(im, fig, ax = None, title=None):
		divider = make_axes_locatable(ax)
		cax = divider.append_axes('right', size='5%', pad=0.1)
		
		cbar = fig.colorbar(im, cax=cax, orientation='vertical')
		cbar.ax.tick_params(labelsize=15)
		if title is not None:
			cbar.ax.set_ylabel(title, rotation=270, fontsize = 18, labelpad=30, y=0.4)

 

	 