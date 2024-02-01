# New POLArized RadIative Simulator NewPOLARIS
** --------------------------- --------------------------- POLARIS code  ---------------------------  ---------------------------  ------------------- 

POLARIS is the 3D radiative transfer code developed by Dr. Stefan Reissl, which is used to: 
    
    . simulate the radiative transfer inside dusty environments, including both dust scattering, absorption, and re-thermal dust emission
    
    . perform synthetic multiwavelength modeling of polarized thermal emission from magnetically aligned dust grains by RAdiative Torques (RATs) and David-Greenstein mechanism.
    
    . simulate self-scattering of thermal dust emission, line radiative transfer, and synchrotron emission.

The link to the original POLARIS code can be found here:
https://github.com/polaris-MCRT/POLARIS





** --------------------------- --------------------------- What is updated in NewPOLARIS  ---------------------------  ---------------------------  -------- 

NewPOLARIS is the extended version from POLARIS including the detailed calculations of the grain alignment by RATs/enhanced Magnetic RAdiative Torques (MRAT) and grain disruption by RAdiative Torque Disruption (RATD). NewPOLARIS sets the platform for connecting grain magnetic properties/internal structure - theory of grain alignment/disruption - synthetic modeling of polarized dust emission - observations of dust polarization.

Properties of aligned dust grains in NewPOLARIS are self-consistently determined based on the grain magnetic properties and conditions of gas density and magnetic fields from input environments. NewPOLARIS will determine:
  
      .  Internal alignment state by Barnett relaxation mechanism

      .  Minimum and maximum alignment size based on the suprathermal rotation and Larmor precession condition

      .  Alignment mechanism (RATs or MRAT) and the corresponding fraction of grains aligning with magnetic fields at high-J attractors  

The Radiative Torques Disruption (RATD) is implemented fully in NewPOLARIS. Given the grain compactness and grain magnetic properties, NewPOLARIS determines:
  
      . Minimum and maximum disruption size
  
      . Fraction of grains being destroyed by RATD and New grain size distribution of enhanced small grains.

      . New dust extinction efficiency, absorption efficiency, and emissivity used for model thermal and polarized dust emission under the change of size distribution by RATD.



** --------------------------- --------------------------- Reference  ---------------------------  ---------------------------  ------------------------ 
POLARIS+ is built based on the theory of grain alignment by RATs (https://academic.oup.com/mnras/article/378/3/910/1052479), MRAT (https://iopscience.iop.org/article/10.3847/0004-637X/831/2/159)), and RATD (https://www.nature.com/articles/s41550-019-0763-6). The first revisit grain alignment theory is done for protostellar cores and disk by Hoang et al.(2022)(https://iopscience.iop.org/article/10.3847/1538-3881/ac9af5). And detailed calculations in POLARIS+ and their applications for understanding the grain alignment and dust polarization in protostellar environments can be found in Giang et al. (2023): https://academic.oup.com/mnras/article/520/3/3788/7044629




 ** --------------------------- --------------------------- How to use it?  ---------------------------  ---------------------------  -------------------- 
To use NewPOLARIS, users need to replace /src and /lib from the original version of POLARIS by /src and /lib of NewPOLARIS, then type:
./compile.sh -u
to compile the new src code.

Please refer to README.cmd for commands used to run new features of grain alignment and disruption.







