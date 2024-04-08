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
# Modelling dust polarization - PM grains
#======================================================================

<task> 1/0
# Detector information: 
     #<detector_dust_polar nr_pixel ="Nx*Ny">  wave_min[m]   wave_max[m]   Nwave   1   angle_from_axis1[degree]   angle_from_axis2[degree]    distance to observer[m]
 
	 <detector_dust_polar nr_pixel ="250*250">	 2000e-6   2000e-6   1    1	   90.00    0.00     3.085678e+18
     <detector_dust_polar nr_pixel ="250*250">   870e-6    870e-6    1    1    90.00    0.00     3.085678e+18
     <detector_dust_polar nr_pixel ="250*250">   450e-6    450e-6    1    1    90.00    0.00     3.085678e+18
 

# Backgorund source information:
 	 <source_dust  nr_photons = "1e6" >
 	 <source_background nr_photons = "1e6" >	 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
 

# Simulation name:
 	 <cmd> CMD_DUST_EMISSION
 

# Directory to the input grid file of previous RAT calculation
 	 <path_grid> "directory/to/POLARIS/example/combine/MRAT/PM/grid_temp.dat"


# Directory to the output folder
 	 <path_out>  "directory/to/POLARIS/example/polarization/MRAT/PM/"
 

# Grain alignment mechanism:
  	 <align> 	 ALIG_RAT
  	 <align> 	 ALIG_INTERNAL
 
 
# Correlation factor between IA and EA:
  	 <f_c> 	0.6


# Fraction of grains at high-J attractor:
     <f_highJ>  0.25

  
 # Account superparamagnetic relaxation effect:
     <change_f_highJ>    1

  
 # Account for inefficient IA due to slow internal relaxation
     <wrong_g_zeta_high_J>   0.15
     <wrong_g_zeta_low_J>    0.05  # If grains at low-J have right IA

     # If grains at low-J with slow internal relaxation have wrong IA, 
     # type, i.e.:
     #<wrong_g_zeta_low_J>    -0.1



 
  # Self-scattering 
 	<rt_scattering> 0
 </task> 
------------------------------------------------------------------------------------------------------------------------------------------------------

 

#======================================================================
# Modelling dust polarization - SPM grains
#======================================================================

<task> 1/0
# Detector information: 
     #<detector_dust_polar nr_pixel ="Nx*Ny">  wave_min[m]   wave_max[m]   Nwave   1   angle_from_axis1[degree]   angle_from_axis2[degree]    distance to observer[m]
 
     <detector_dust_polar nr_pixel ="250*250">   2000e-6   2000e-6   1    1    90.00    0.00     3.085678e+18
     <detector_dust_polar nr_pixel ="250*250">   870e-6    870e-6    1    1    90.00    0.00     3.085678e+18
     <detector_dust_polar nr_pixel ="250*250">   450e-6    450e-6    1    1    90.00    0.00     3.085678e+18
 

# Backgorund source information:
     <source_dust  nr_photons = "1e6" >
     <source_background nr_photons = "1e6" >     0.0 0.0 0.0 0.0 0.0 0.0 0.0 
 

# Simulation name:
     <cmd> CMD_DUST_EMISSION
 

# Directory to the input grid file of previous RAT calculation
      <path_grid> "directory/to/POLARIS/example/combine/MRAT/SPM_1e4/grid_temp.dat"


# Directory to the output folder
      <path_out>  "directory/to/POLARIS/example/polarization/MRAT/SPM_1e4/"
 

# Grain alignment mechanism:
     <align>     ALIG_RAT
     <align>     ALIG_INTERNAL
 
 
# Correlation factor between IA and EA:
     <f_c>  0.6


# Fraction of grains at high-J attractor:
     <f_highJ>  0.25

  
 # Account superparamagnetic relaxation effect:
     <change_f_highJ>    1

  
 # Account inefficient IA effect
     <wrong_g_zeta_high_J>   0.15
     <wrong_g_zeta_low_J>    0.05

     # If grains at low-J with slow internal relaxation have wrong IA, 
     # type, i.e.:
     #<wrong_g_zeta_low_J>    -0.1
     
 
  # Self-scattering 
    <rt_scattering> 0
</task> 
------------------------------------------------------------------------------------------------------------------------------------------------------

 

 