# EXTPAR

EXTPAR (External Parameters for Numerical Weather Prediction and Climate Application) is an official software of the COSMO Consortium (www.cosmo-model.org).  It is used to prepare the external parameter data files that are used as input for the COSMO model, and additionally now the ICON model.    

Further information can be found in the [manual](doc/user_and_implementation_manual.pdf)

For compilation instructions see: [README.compile_run](doc/README.compile_run)

# Release notes

### 5.0
This release represents a merge of the Extpar official version 4.0 code with the DWD-Extpar version 2.10.  

#### Build Mechanism

* Added Options file for compiling on LCE with Intel compiler 
* Added Options files for compiling on Mistral with GCC, NAG, and Intel compilers
* Updated Options file for compiling on o3 at ETH with PGI compiler
* Update bin/gen_info.sh for compatibility with git

#### Run scripts

* Minor adaptations to run scripts to be compatible with new version
* Addition of MPI ICON run script

#### Code changes

#### Albedo calculation:
* Added bin/cdo2alb-buffer.py to replace slow albedo calculation for ICON

#### Topography calculation: 
* Added namelist parameter: subtract mean_slope
* Move mo_agg_topo to mo_agg_topo_icon and mo_agg_topo_cosmo

#### Soil calculation: 
* Additional HWSD calculation for deep soil

#### NDVI calculation:
* Added bin/cdo2ndvi-buffer.py to replace slow NDVI calculation for ICON

#### Climatological 2M temperature calculation:
* Changed namelist parameter it_cl_type to raw_data_t_id

#### Flake calculation:
* Added namelist parameter: lflake_correction

#### Consistency check: 
* [Changes results] Included lower limit for roughness length
* ERA-I SST and T2M temperature for ICON model

#### Grib output NOT supported

#### GME model NOT supported


