# POLARIS (POLArized RadIation Simulator) + MAgnetic Radiative Torque Alignment


In this version, we incorporate the grain magnetic properties and the detailed internal and external alignment into POLARIS. The new in this version is as follows:

1. We separate CMD_RAT from CMD_TEMP.

Previous, user can model the grain alignment via <cmd> CMD_TEMP_RAT or <cmd> CMD_RAT. But we found that even with CMD_RAT, POLARIS still needs to calculate the mcrt and dust temperature before
modelling RAT alignment. It is inconvenient when we want to model the grain alignment with different grain magnetic properties. Therefore, we create the option to run CMD_RAT directly from the output file from CMD_TEMP. 

In detail, with CMD_TEMP, we write the radiation field into the grid file (grid_temp.dat), so with CMD_RAT, we can read information of the radiation field from grid_temp.dat to calculate Raditiative torque and align. 






2.  Users can define the grain type (PM or SPM grains) and the configuration of iron inside dust grains (for MRAT alignment only):

	- For PM grains, iron will be distributed diffusely inside grains with the iron fraction f_p is input by hand via the command:

		<iron_fraction> 0.1 [this value can be changed)
		
	- For SPM grains, iron atoms will be under the cluster form described by the number of iron/cluster and the volume filling factor of iron clusters. These information are inputted with:

		<number_cluster> 1e2

		<volume_filling_factor> 0.1





3. We separate the alignment modelling with RATs and MRAT alignment.

	- For RAT alignment, type <align> ALIG_RAT. In this mode:
	
		+ We only calculate the minimum alignment size [a_align] defined by the suprathermal rotation by RAdiative Torques. 
		+ The maximum alignment size [amax_JB_Lar] is set to be equal the chosen maximum grain size. 
	
	
	- For MRAT alignmentm type <align> ALIG_MRAT. In this mode:

		+ We model the internal and external alignment self-consistently with the input grain magnetic properties. To be more detail:
		
		+ For the internal alignment, we first calculate the Barnett relaxation timescale based on the input information of the grain magnetic properties. 
		  Then based on the condition tau_bar <= tau_gas, we the minimum and maximum size for grains having fast internal relaxation at low and high-J attractor 
			  . parameter [amin_aJ_lowJ] - [amax_aJ_lowJ] and [amin_aJ_highJ] - [amax_aJ_highJ]. 
			 
			 
		+ For the external alignment:
		
			- We first calculate the minimum alignment size [align] based on RAT criteria.
			
			- We then calculate the Larmor precession timescale and base on the condition tau_Lar <= tau_gas/10 (Yang et al. (2021)) to determine the maximum alignment size [amax_JB_Lar].
			
			- For the value of f_highJ, we calculate the magnetic relaxation timescale. We base on the condition of the magnetic relaxation parameter delta_m (see Hoang+2016, Giang+2022) to determine:
				- the minimum and maximum size for 50% of grains will be aligned with B at high-J attractor [aminJB_DG_0.5] and [amaxJB_DG_0.5]
				- the minimum and maximum size for 100% of grains having perfect alignment [aminJB_DG_1] and [amaxJB_DG_1]
			
				* Indeed, currently we assume f_highJ = 0.5 for grains having 1 <= delta_m <= 10. This value can be changed depending on the purpose of users. 
				To do it, change the line 3156 in file Dust.cpp, i.e., from f_highJ = 0.5; to f_highJ = [value you want]. 
			 





4. For the calculation of the Rayleigh reduction factor R, we dont need to add the command <align> ALIG_INTERNAL to calculate the internal alignemnt degree of grains.

	- For <align> ALIG_RAT, users only need to provide the value for <f_highJ> to calculate R
	
	- For <align> ALIG_MRAT:
	
		+ To activate the dependence of f_highJ with grain size as described in line [47-48], type <change_f_highJ> 1 [in default, we do not turn on this option, or <change_f_highJ> 0]
	
		+ For the internal alignment degree of grains having slow internal relaxation, type the command <wrong_g_zeta_low_J> for grains at low-J, and 
		  <wrong_g_zeta_high_J> for grains at high-J
			 

 

## Installation

Open a terminal/console and move into the POLARIS directory:
```bash
cd /YOUR/POLARIS/PATH/
```

To install POLARIS, run the installation script:
```bash
./compile.sh -f
``` 


We attach the example file how to active MRAT calculation in the /example (for both CMD_TEMP, CMD_RAT, and CMD_DUST_EMISSION)
