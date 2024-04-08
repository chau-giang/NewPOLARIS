#======================================================================
# Header
#======================================================================
<common> 
# Dust information
  <dust_component>  "/directory/to/POLARIS/input/dust/silicate_oblate.dat" "plaw" 0.625 3800 5e-09 1e-4 -3.5
  <dust_component>  "/directory/to/POLARIS/input/dust/graphite_oblate.dat" "plaw" 0.375 2250 5e-09 1e-4 -3.5

  <mass_fraction>       0.01

    <axis1> 1   0   0  # x axis
    <axis2> 0   1   0  # y axis

    <write_inp_midplanes>   256   # Number of pixel in input fits file
    <write_out_midplanes>   256   # Number of pixel in output fits file
  <write_radiation_field> 1   # To write radiation field into the fits file. [0]: dont write it.

    <mu>            2.0               # Average atomic number

    <phase_function>    PH_ISO

    <nr_threads>    -1              # to use all available cores in the computer
</common>

#======================================================================
# Modelling dust polarization
#======================================================================

<task> 1
 # Detector information: 
     #<detector_dust_polar nr_pixel ="Nx*Ny">  wave_min[m]   wave_max[m]   Nwave   1   angle_from_axis1[degree]   angle_from_axis2[degree]    distance to observer[m]
 
	 <detector_dust_polar nr_pixel ="250*250">	 2000e-6   2000e-6   1    1	   90.00    0.00     3.085678e+18
     <detector_dust_polar nr_pixel ="250*250">   870e-6    870e-6    1    1    90.00    0.00     3.085678e+18
     <detector_dust_polar nr_pixel ="250*250">   450e-6    450e-6    1    1    90.00    0.00     3.085678e+18
 

 # Backgorund source information:
     #<source_dust  nr_photons = "Nph" > 
 	 <source_dust  nr_photons = "1e6" >
 	 <source_background nr_photons = "1e6" >	 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
 


 # Simulation name:
 	 <cmd> CMD_DUST_EMISSION
 


 # Directory to the input grid file of previous RAT calculation
 	 <path_grid> "directory/to/POLARIS/example/combine/RATs/grid_temp.dat"


# Directory to the output folder
 	 <path_out>  "directory/to/POLARIS/example/polarization/RATs/"
 

 # Grain alignment mechanism:
  	 <align> 	 ALIG_RAT
  	 <align> 	 ALIG_INTERNAL
 
 
  # Correlation factor between IA and EA:
  	 <f_c> 	0.6



 # Fraction of grains at high-J attractor:
  	 <f_highJ> 	1.00
 

 
  # Self-scattering 
 	<rt_scattering> 0
 </task> 
------------------------------------------------------------------------------------------------------------------------------------------------------

 


 