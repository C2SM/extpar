# Release notes
 ## 5.3
 This is an intermediate release that reduces code complexity for topo_to_buffer.exe, enhances the testing for INTEL compiler and further cleans the code

 * Merge sgsl_to_buffer into topo_to_buffer
   - New namelist &oro_runcontrol in INPUT_ORO containing lcompute_sgsl
   - For users of former sgsl_to_buffer.exe, namelist &sgsl_io_extpar now moved to INPUT_ORO, containing the new parameter lpreproc_oro
   - The functionality is kept by default for all newly introduced namelist parameters, so for the same workflows as before only change lcompute_sgsl

 * Testsuite
   - Additional check for compiler warning of GCC,INTEL and NAG
   - Unify runscripts for COSMO and ICON
   - Slightly different domain for COSMO1, reducing the required number of ASTER tiles to only 1
   - Add references for INTEL in a seperate directory in data

 * Cleanup
   - Initialize logicals in extpar_consistency_check properly to prevent bugs
   - Remove hardcoded filename in emiss_to_buffer
   - Finalize logging and coding standard as described in developers guide

 ## 5.2.1
 This is a minor release containing a bug fix and a small feature addition.
 * Bug fix for ICON/COSMO file- and variable name mismatch in topography calculation
 * Add Extpar version number (pulled from git release number) to output NetCDF file
 
 ## 5.2
 This is an intermediate release introducing extpar_emiss_to_buffer, an improved logging, enhanced error checking during I/O and a lot of clean-up and formatting

 * New Extpar executable emiss_to_buffer
   - Aggregates CAMEL emissivity data to the target grid
   - Two raw datasets available (full range and only long-wave radiation)

 * Consistent logger for all Extpar executables
   - Three levels of messages: info, warning and error
   - Each Extpar executable write to its own logfile

 * Clean-up and formatting of all src-files
   - Remove all unused variables and USE-statements
   - Remove all unused dummy arguments in subroutines
   - Implement formatting according the coding-guidelines for Extpar
 
 * Make all precisions consistent
   - Remove i8 from Extpar, instead make all INTEGER(KIND=i4)
   - Change all REAL to REAL(KIND=wp), wp is defined in mo_kind
 
 * Small changes in some fields due to fix of implicit type conversion during runtime
   - Z0, max difference ~10^(-7)
   - DEPTH_LK, max difference ~10^(-6)

 * Enhanced error checking during I/O
   - All namelist I/O checked, abort of Extpar in case of incorrect (typos, wrong variables, etc,) namelists


 ## 5.1.2
 This is a minor release containing a few bug fixes.
 * Fix build environment on Kesch
 * Add missing definition of skinc_lu meta data when ECOCLIMAP dataset is chosen. 

 ## 5.1.1
 This is a minor release containing a few bug fixes.
 * Fix read of l_use_glcc landuse calculation for COSMO runs.
 * Fix unitialized logical flag to trigger scale separation in topography calculation.
 * Reactivate all cosmo tests from testsuite on Kesch.
 
 ## 5.1
 This is an intermediate release containing some bug fixes and some minor developments.  

 * Changes to Jenkins and the automated testing:  
   - Fix Mistral setup so that code runs on compute nodes instead of login nodes
   - Fix Mistral setup so that intel compiler can be tested.  Note that only run success
     is checked for the intel compiler; the results are not yet tested.  
   - Fix NAG compiler setup so that only compilation, not testing is done, because
     testing is too time consuming.  

 * Bug fix for iaot_type = 4 (MACv2 aerosols).  The code had not been correctly imported from version 4.0.
 
 * Contributions from DWD including:
   - DWD versions of the python and shell replacement scripts for ndvi, albedo, and cru
   - DWD bug fix for the albedo calculation
   - DWD bug fix for incorrect glacier points

 * New output variable skin conductivity (SKC) developed by Jan-Peter Schulz.  Skin conductivity is calculated from the landuse data.  

 ## 5.0.4
 This is a minor release containing a few bug fixes.

 * Bug fix for problems when soil_type=3 is used.  The code had not been correctly imported from version 4.0  
 
 * Bug fix adding missing NetCDF get_varid call when more than one GLOBCOVER tile is used.  

 ## 5.0.3
 This is a minor release containing a bug fix.  

 * Bug fix for incorrect global attributes in output NetCDF file.  

 ## 5.0.2
 This is a minor release containing a bug fix.

 * Bug fix for incorrect subroutine argument usage in mo_agg_isa.f90.  

 ## 5.0.1
 This is a minor release fixing a few bugs and some documentation.

 * Bug fix to remove too many characters in write statement (mo_logging.f90)
 * Bug fix for unitialized variable ntiles_globcover (mo_landuse_routines.f90, extpar_consistency_check.f90, and extpar_landuse_to_buffer.f90)
 * Remove unnecessary libraries in Options.daint
 * Fix typos in README.compile_run

 ## 5.0
 This release represents a merge of the Extpar official version 4.0 code with the DWD-Extpar version 2.10. 

 #### Build Mechanism

 * Added Options file for compiling on LCE with Intel compiler
 * Added Options files for compiling on Mistral with GCC, NAG, and Intel compilers
 * Added Options files for compiling on Mac OS with the GCC compiler
 * Updated Options file for compiling on o3 at ETH with PGI compiler
 * Update bin/gen_info.sh for compatibility with git

 #### Run scripts

 * Minor adaptations to run scripts to be compatible with new version
 * Addition of MPI ICON run script

 #### Testing

 * Addition of cosmo-dwd and icon tests to testsuite.
 * Update of Jenkins build and test scripts.
 * Add testsuite run scripts for mpi and dwd on mistral.
 * Update of testsuite references.

 ### Code changes

 #### Albedo calculation:
 * Added bin/cdo2alb-buffer.py and extpar_alb_to_buffer.sh to replace slow and incorrect albedo calculation for high resolution ICON model grids

 #### Topography calculation:
 * Added namelist parameter: lsubtract mean_slope
 * Move mo_agg_topo to mo_agg_topo_icon and mo_agg_topo_cosmo
 * Bug fix in mo_topo_sso- changes results of sso_sigma slightly
 
 #### Soil calculation:
 * Additional HWSD calculation for deep soil
 
 #### NDVI calculation:
 * Added bin/cdo2ndvi-buffer.py and extpar_ndvi_to_buffer.sh to replace slow and incorrect NDVI calculation for high resolution ICON model grids
 
 #### Climatological 2M temperature calculation:
 * Added bin/cdo2t_cl-buffer.py and extpar_cru_to_buffer.sh to replace slow and incorrect TCLIM calculation for high resolution ICON model grids
 * Added namelist parameter: ltcl_merge
 
 #### Flake calculation:
 * Added namelist parameter: lflake_correction.
 
 #### Consistency check:
 * [Changes results] Included lower limit for roughness length
 * ERA-I SST and T2M temperature for ICON model

 #### Grib output NOT supported

 #### GME model NOT supported

 ### Changes in Results
  Due to the large amount of changes in the code in this release, there are many differences in the resulting external parameter fields generated by the release 5.0 code compared to the fields generated by older Extpar codes.  The only code change in this release that deliberately changed the results was the addition of a lower limit for roughness length, which is 1e-6.  Otherwise, any changes in results that can be seen came directly from bug fixes to the code, and as such most of them are small and are not expected to change results in the COSMO or ICON model runs.  Some of these changes in results are examined in the next two sections.  

 #### Extpar version 4.0 to Extpar version 5.0
  The technical testsuite in Extpar was used to compare the external parameter fields from Version 4.0 and Version 5.0 for three different MeteoSwiss operational setups for COSMO and a climate setup for COSMO-CLM.  

 ##### COSMO 7, globe topography input
  For the COSMO7 MCH setup using the globe topography data set, changes in the results smaller than 1e-6 can be seen for several variables, including the aerosol variables, albedo variables, HORIZON, and SSO_SIGMA.  Larger changes on the order of 3 degrees can be seen in the SSO_THETA variable;  these are due to a bug fix in this release, and are expected.  Finally, roughness length is different as well due to the introduction of the lower limit value of 1e-6.  

 ##### COSMO 7, aster topography input
  For the COSMO7 MCH setup using  the aster topography data set, the changes in the results are less than 5e-7, and occur in the aerosol, HORIZON, and SSO_SIGMA variables.  

 ##### COSMO 1, aster topography input
  For the COSMO1 MCH setup using the aster topography data set, the changes in the results are less than 4e-5, and occur in the aerosol, HORIZON, SKYVIEW, and T_CL variables.  

 ##### COSMO-CM climate setup
  For the COSMO-CLM climate setup using the globe topography data set, changes in results smaller than 3e-8 can be seen in the SSO_SIGMA and SSO_STDH variables.  Larger changes on the order of 3 degrees can be seen in the SSO_THETA variable;  these are due to a bug fix in this release, and are expected.  Due to another bug fix, the ALB_SAT and ALB_DRY variables have changed results on the order of .2.  Finally, roughness length is different as well due to the introduction of the lower limit value of 1e-6.  

 #### DWD Extpar version 2.10 to Extpar version 5.0
  Comparisons of the COSMO D2 setup used operationally by DWD were carried out to compare the current operational Extpar code (DWD version 2.10) with the new release 5.0.  This comparison showed no significant differences in the generated external parameter fields.  Get more details of this comparison [here](doc/EXTPAR_CD2_comparison.pdf)
