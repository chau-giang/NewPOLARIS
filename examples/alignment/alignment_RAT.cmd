#======================================================================
# Header
#======================================================================
<common> 
# Dust information
  <dust_component>  "/directory/to/POLARIS/input/dust/silicate_oblate.dat" "plaw" 0.625 3800 5e-09 1e-4 -3.5
  <dust_component>  "/directory/to/POLARIS/input/dust/graphite_oblate.dat" "plaw" 0.375 2250 5e-09 1e-4 -3.5

  <mass_fraction>		0.01

	<axis1>	1	0	0
	<axis2>	0	1	0

	<write_inp_midplanes>	256   # Number of pixel in input fits file
	<write_out_midplanes>	256   # Number of pixel in output fits file
  <write_radiation_field> 1   # To write radiation field into the fits file. [0]: dont write it.

	<mu>			2.0               # Average atomic number

	<phase_function>	PH_ISO

	<nr_threads>	-1          	# to use all available cores in the computer
</common>

 
#======================================================================
# Calculation of grain alignment by RATs
#======================================================================
<task> 1
    # Radiation source
  	   <source_star nr_photons = "5e4">	0	0	0   2  6000
  	   <source_isrf nr_photons = "5e4"> 1 1
  

    # Calculation of grain alignment distribution
  	<cmd>			CMD_RAT 
  

    # Directory to the input grid file of previous MCRT calculation
  	 <path_grid>		"directory/to/POLARIS/example/MCRT/grid_temp.dat"

    # Directory to the output folder
	   <path_out>		"directory/to/POLARIS/example/alignment/RATs/"
  	 
    # Setting calculation mode
     <dust_offset>		 1

  # Magnetic properties of grains:
   <align>  ALIG_RAT
      
  # Unit conversion: # to convert to SI unit.
  	<conv_dens>		1.0
  	<conv_len>		1.0
  	<conv_mag>		1.0
  	<conv_vel>		1.0
</task> 
    
    
    
  