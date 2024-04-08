
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

 

class fits_file_pol_map:
	def __init__(self, filename):
		self.filename = filename

	def read_fits(self):
		"""
		Args:
			filename: directory to fits file
	 
		Returns:
			dictionary of fits file
		"""

		hdulist = fits.open(self.filename)
		hdu = hdulist[0]

		dat = dict()
		dat_ = hdu.data
		head_ = hdu.header

		#------------------------------------------------------------------------
		dat['boundary'] = dict()
		dat['boundary']['Nx'] = head_['NAXIS1']
		dat['boundary']['Ny'] = head_['NAXIS2']
 
		dat['boundary']['dx'] = head_['CDELT1B']
		dat['boundary']['dy'] = head_['CDELT2B']
		dat['boundary']['x'] = head_['CRVAL1B'] + np.arange(dat['boundary']['Nx']) * dat['boundary']['dx']
		dat['boundary']['y'] = head_['CRVAL2B'] + np.arange(dat['boundary']['Ny']) * dat['boundary']['dy']
		dat['boundary']['xmin'] = head_['CRVAL1B'] - 0.5 * dat['boundary']['dx']
		dat['boundary']['ymin'] = head_['CRVAL2B'] - 0.5 * dat['boundary']['dy']
		dat['boundary']['xmax'] = dat['boundary']['x'].max() + 0.5 * dat['boundary']['dx']
		dat['boundary']['ymax'] = dat['boundary']['y'].max() + 0.5 * dat['boundary']['dy']
		dat['boundary']['angle_axis1'] = head_['RANGLE1']
		dat['boundary']['angle_axis2'] = head_['RANGLE2']

		#------------------------------------------------------------------------
		dat['data'] = dict()
		dat['data']['nr_wave'] = head_['NAXIS3']

		if dat['data']['nr_wave'] == 1:
			dat['data']['wavelength'] = head_['HIERARCH WAVELENGTH1']
		else:
			wave = np.zeros(dat['data']['nr_wave'])
			for i in range (0, dat['data']['nr_wave']):
				wave[i] =  head_['HIERARCH WAVELENGTH'+str(i+1)]
			dat['data']['wavelength'] = wave

		dat['data']['Stoke_I'] = dat_[0, :, :, :]
		dat['data']['Stoke_Q'] = dat_[1, :, :, :]
		dat['data']['Stoke_U'] = dat_[2, :, :, :]
		dat['data']['Stoke_V'] = dat_[3, :, :, :]
		dat['data']['P'] = np.sqrt(pow(dat['data']['Stoke_Q'], 2) + pow(dat['data']['Stoke_U'], 2)) / dat['data']['Stoke_I']
		dat['data']['PA_rad'] = np.arctan2(dat['data']['Stoke_U'], dat['data']['Stoke_Q'])/2
		dat['data']['PA_degree'] = dat['data']['PA_rad']*180/math.pi
 
		dat['data']['tau'] = dat_[4, :, :, :]
		dat['data']['NH'] = dat_[5, :, :, :]*1e-4 #cm-2
		return dat


	def plot_polarization_map(self, fig, ax=None, model_name = '', step_reduce = 15, font_size = 25, units = 'width', scale_units = 'width', 
		color = 'black', width = 0.006, headwidth = 1, headlength = 1,	headaxislength=6,	pivot = 'middle', linewidth = 1, title = True, xlabel=True, 
		ylabel=True, titcolor=True, plot_B_vector=False, plot_intensity=False, only_vector=False):
		"""  
		Args:
			Direction to filename
			step_reduce: size of reduced grid to visualize the polarization vector
			plot_B_vector: [False]: plot polarization (E) vector, [True]: plot polarization (B) vector
			plot_intensity: [False]: color code shows polarization degree. [True]: color code shows intensity
			only_vector: [False]: vector length present polarization degree. [True]: same length, only show polarization vector direction.

		Return:
			Polarization map with the polarization vector
		"""

		if ax is None:
			plt.gca()

		dictionary = self.read_fits()
		I = dictionary['data']['Stoke_I'][0, :, :]
		P = dictionary['data']['P'][0, :, :]*100

		# Set nan pixel to min(P)
		index_nan = np.where(np.isnan(P))
		nr_pixel_nan = len(index_nan[0])
		P[index_nan[0], index_nan[1]] = np.nanmin(P)
		
		xmin = dictionary['boundary']['xmin']
		xmax = dictionary['boundary']['xmax']
		ymin = dictionary['boundary']['ymin']
		ymax = dictionary['boundary']['ymax']
		Nx = dictionary['boundary']['Nx']
		Ny = dictionary['boundary']['Ny']
		boundary = [xmin, xmax, ymin, ymax]
		if plot_intensity: # If plot intensity map
			im = ax.imshow(I, origin = 'lower', cmap = 'viridis', norm = matplotlib.colors.LogNorm(), extent = boundary)	
		else: # Plot polarization degree map
			im = ax.imshow(P, origin = 'lower', cmap = 'viridis',
				vmin = np.nanmin(P), vmax = np.nanmax(P), extent = boundary)

		if titcolor == True:
			if plot_intensity: 
				configure_colorbar(im, fig, ax, title = 'Intensity')
			else:
				configure_colorbar(im, fig, ax, title = r'$\sf{\rm P(\%)}$')
		else:
			titcolor = True
			configure_colorbar(im, fig, ax)
 
		avg_angle = dictionary['data']['PA_rad'][0, :, :]
		if plot_B_vector: #if plot polarization (B) vector:
			avg_angle = avg_angle + np.pi/2

		if only_vector: # if only plot polarization vector direction
			P_x = np.cos(avg_angle)
			P_y = np.sin(avg_angle)
		else:
			P_x = P * np.cos(avg_angle)
			P_y = P * np.sin(avg_angle)

		x = np.linspace(xmin, xmax, Nx)
		y = np.linspace(ymin, ymax, Ny)
		X, Y = np.meshgrid(x, y)


		# To determine the color of polarization vector. If almost pixel has P > <P> -> pol vec should be black
		#												 Otherwise, pol vec will be white
		nr_pixel_notnan = dictionary['boundary']['Nx']*dictionary['boundary']['Ny'] - nr_pixel_nan

		average_P = (np.nanmin(P) + np.nanmax(P))/2
		index_above_aveP = np.where(P > average_P)[0]
		if len(index_above_aveP) > 0.5*(nr_pixel_notnan//2):
			colors = 'white'
		else:
			#colors = 'black'
			colors = 'white'


		ax.quiver(X[::step_reduce, ::step_reduce], Y[::step_reduce, ::step_reduce], P_x[::step_reduce, ::step_reduce],
		 		  P_y[::step_reduce, ::step_reduce], color = colors, units = 'width', scale_units = 'width', width = width,
		 		headwidth = headwidth,	headlength = headlength, headaxislength=headaxislength,	pivot = 'middle',
		 		linewidth = linewidth)


		# Configure figure
		configure_figure(ax)

		# Label for x axis
		if xlabel == True:
			ax.set_xlabel('x [au]', fontdict={'fontsize': 16})

		# Label for y axis
		if ylabel == True:	
			ax.set_ylabel('z [au]', fontdict={'fontsize': 16})
		
		# Title of Figure
		wavelength = dictionary['data']['wavelength']*1e6
		if title == True:
			# Model name and wavelength
			props = dict(boxstyle='round', facecolor='white', alpha=0.9)
			# place a text box in upper left in axes coords
			ax.text(0.03, 0.98, ' %s %d %s' % (model_name, wavelength, r'$\sf{\mu m}$'), transform=ax.transAxes, fontsize=13.5,
					color = 'black', verticalalignment='top', bbox=props)
 
		return(ax)
 

	def index_radius(self, dictionary, distance, nr_pixel_beam_size):  
		"""
		Args: 
			dictionary: dictionary of fits file
			distance: radius from the center region
			nr_pixel_beam_size: number of pixel in the beam size 

		Returns:
			coor_x: position of cell in the observed region on x direction
			coor_y: position of cell in the observed region on y direction
		"""
		nr_pixel = dictionary['boundary']['Nx']
		x = np.arange(0, nr_pixel, 1)   
		X, Y = np.meshgrid(x, x)
		R = np.sqrt((X - nr_pixel//2)**2 + (Y - nr_pixel//2)**2)

		if ((distance != 0) and (distance - nr_pixel_beam_size > 0)): 
			index = np.where(R < distance)		# remove grid inside the distance
			coor_x = index[0]
			coor_y = index[1]
			R[coor_x, coor_y] = np.nan

		#if ((distance != (nr_pixel//2)-1) and (distance + nr_pixel_beam_size < nr_pixel)):
		index = np.where(R > distance+nr_pixel_beam_size)	#remove grid outside the distance
		coor_x = index[0]
		coor_y = index[1]
		R[coor_x, coor_y] = np.nan	

		index = np.where(np.isnan(R) == False)
		coor_x = index[0]
		coor_y = index[1]
		return coor_x, coor_y


	def P_I_matrix(self, dictionary, nr_pixel_beam_size = 1):           
		"""
		This function is to get the 
		Args:
			dictionary: dictionary of fits file
			nr_pixel_beam_size: number of pixel within the beam size, default is 1 pixel

		Returns:
			2D matrix of [nr_pixel_in_outer_boundary x nr_pixel_radial_direction]
			of I, P, tau, and NH
			number of value (# 0) on each column
		"""

		max_grid = 0
		nr_pixel = dictionary['boundary']['Nx']//2
		for i in range(0, nr_pixel):
			coor_x, coor_y = self.index_radius(dictionary, i, nr_pixel_beam_size)
			if (len(coor_x) > max_grid):
				max_grid = len(coor_x)

		I_matrix = np.zeros([max_grid, nr_pixel])   
		P_matrix = I_matrix.copy()                  		                                           
		tau_matrix = I_matrix.copy()
		NH_matrix = I_matrix.copy()
		nr_point = np.zeros(nr_pixel)

		for i in range(0, nr_pixel):	
			coor_x, coor_y = self.index_radius(dictionary, i,  nr_pixel_beam_size)

			I = dictionary['data']['Stoke_I'][0, coor_x, coor_y]
			P = dictionary['data']['P'][0, coor_x, coor_y]
			tau = dictionary['data']['tau'][0, coor_x, coor_y]
			NH = dictionary['data']['NH'][0, coor_x, coor_y]*1e-4 # cm-2

			cout = 0
			for j in range(0, len(P)):
				I_matrix[j,i] = I[j]
				P_matrix[j,i] = P[j]				
				tau_matrix[j,i] = tau[j]
				NH_matrix[j,i] = NH[j]
				cout += 1
			nr_point[i] = cout
                    
		return I_matrix, P_matrix, tau_matrix, NH_matrix, nr_point 


	def mean_value(self, I_matrix, P_matrix, tau_matrix, NH_matrix, nr_point):
		"""
		Args: 
			2D matrix of [nr_pixel_in_outer_boundary x nr_pixel_radial_direction] of
			I, P, tau, NH
			nr of point (value # 0) in each column of the matrix

		Returns:
			array of the mean value over the circle at different position from the center 
			to the outer boundary
		"""
		I_mean = np.zeros(len(I_matrix[0, :]))
		P_mean = I_mean.copy()
		tau_mean = I_mean.copy()
		NH_mean = I_mean.copy()
		for i in range(0, len(I_mean)):
			I_mean[i] = sum(I_matrix[:, i])/nr_point[i]
			P_mean[i] = sum(P_matrix[:, i])/nr_point[i]
			tau_mean[i] = sum(tau_matrix[:, i])/nr_point[i]
			NH_mean[i] = sum(NH_matrix[:, i])/nr_point[i]
		return I_mean, P_mean, tau_mean, NH_mean		
        

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


    