#======================================================================
# Header
#======================================================================
<common> 
# Dust information
  #<dust_component>   "/directory/to/POLARIS/input/dust/silicate_oblate.dat"  "plaw"  fsil   rho_sil   amin  amax  eta
  #<dust_component>   "/directory/to/POLARIS/input/dust/graphite_oblate.dat"  "plaw"  fcarc  rho_carc  amin  amax  eta

	<dust_component>	"/directory/to/POLARIS/input/dust/silicate_oblate.dat" "plaw" 0.625 3800 5e-09 1e-4 -3.5
	<dust_component>	"/directory/to/POLARIS/input/dust/graphite_oblate.dat" "plaw" 0.375 2250 5e-09 1e-4 -3.5

  <mass_fraction>		0.01

	<axis1>	1	0	0
	<axis2>	0	1	0

	<write_inp_midplanes>	256
	<write_out_midplanes>	256
	<plot_inp_midplanes>	0
	<plot_out_midplanes>	0
  <write_radiation_field> 1 # to save radiation field into output fits file. [0]: dont save.

	<mu>			2.0 # Average atomic number

	<phase_function>	PH_ISO

	<nr_threads>		-1 # to use all available cores in the computer
</common>

 
#======================================================================
# MCRT calculation
#======================================================================
<task> 1
    # Stellar radiation source
      <source_star nr_photons ="Nph" >    x_star y_star z_star R_star T_star
  	  <source_star nr_photons = "5e4">	0	0	0   2  6000

      <source_isrf nr_photons ="Nph" >    G0 R0
  	  <source_isrf nr_photons = "5e4"> 1 1
  

    # Calculation of radiation field and dust/gas temperature
  	  <cmd>			CMD_TEMP  
  


    # Directory to the input grid file
  	  <path_grid>		"directory/to/POLARIS/example/combine/grid.dat"
	  

    # Directory to the output folder
      <path_out>		"directory/to/POLARIS/example/combine/RATs/"



    # Setting calculation mode
      <dust_offset> 1       #To use the default dust temperature in the input grid.dat file. [0]: dont use. 
      <radiation_field> 1   #To save the radiation file into the grid file. [0]: dont use.
      <full_dust_temp> 1    #To calculate and save grain temperature of the full grain size distribution. 
                                [0]: only calculate the average dust temperature over grain size distribution. 
      <sub_dust> 1          #To consider the sublimation of dust grains. [0]: dont use.
      <adj_tgas>    1       #Dust to gas ratio

  
    # Magnetic properties of grains:
      <align>  ALIG_RAT

    # Unit conversion: # to convert to SI unit.
      <conv_dens>   1
      <conv_len>    1
      <conv_mag>    1
</task> 
 
 
    
    
    
   