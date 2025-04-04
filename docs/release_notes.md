# Release notes

For the latest release notes, please see our [GitHub webpage :material-open-in-new:](https://github.com/C2SM-RCM/extpar/releases){:target="_blank"}.

## 5.14

* New buffer script `extpar_cdnc_to_buffer` for ICON
  - Cloud droplet number concentration from the cloud optical depth and effective radius
  - Input data:
    - modis_cdnc_climatology_Q06.nc

## 5.13

* New buffer script `extpar_edgar_to_buffer` for ICON
  - Global emission data for black carbon, organic carbon and sulfur dioxide
  - Input data:
     - EDGARv6.1_BC_2018_TOTALS.0.1x0.1.nc
     - EDGARv6.1_OC_2018_TOTALS.0.1x0.1.nc
     - EDGARv6.1_SO2_2018_TOTALS.0.1x0.1.nc
* Fix albodo for glacier points
   - Changes the following fields for ICON setups
      - ALB
      - ALNID
      - ALUVD
      - FOR_D
      - FOR_E
      - LAI_MX
      - PLCOV_MX
      - ROOTDP
      - SOILTYP
      - ROOTDP
      - URBAN

## 5.12

* New landuse data set Ecoclimap Second Generation
  - Namelist switch `i_landuse_data=6` in `INPUT_LU` 
  - Set `l_terra_urb=.true.` in `INPUT_LU` to process additional urban fields
  - Corresponding input-data set is `ECOCLIMAP_SG.nc`
* Assign more landuse classes (previously "undefined") for Corine
  - Changes the following fields
     - SOILTYP
     - FR_LAND
     - PLCOV_MX
     - LAI_MX
     - RSMIN
     - URBAN
     - FOR_D
     - FOR_E
     - SKC
     - EMIS_RAD
     - ROOTDP
     - Z0
     - NDVI_MAX
     - FR_LAKE
     - DEPTH_LK
     - AHF
     - ISA
     - NDVI
     - ALNID
     - ALUVD
     - NDVI_MRAT
     - LU_CLASS_FRACTION
     - ALB
* New buffer script `extpar_hwsdART_to_buffer`
  - Standalone executable, output not processed by consistency_check
  - Namelist `INPUT_hwsdART` defines key parameters
  - Input data: HWSD0_USDA.nc
  - Problems with git-lfs prohibit the data to be provided on Levante

## 5.11

* Remove ecoclimap (ilanduse=4) from extpar_landuse_to_buffer
* Introduction of GitHub actions
  - Format code each time a .py-file is pushed to repo
  - Test if docs (.tex-file) is still buildable each time it is changed
  - Attach docs (.pdf) as artifact to each git tag
* Remove misleading code for extpar_topo_to_buffer
  - Orographic smoothing for Icon aborts Extpar
  - Remove dead code parts in mo_agg_topo_icon
* Installation of Extpar as Python package
  - Possibility to install all relevant executables and Python scripts as package
  - Adapt import statements in all Python files to work as a package as well
  - Detailed instructions how to install in [README](README.md)
* [Wrapper script](python/WrapExtpar.py) for Extpar
  - Only key switches like _itope_type_ or _ialb_type_ need to be specified
  - Automatic generation of namelists
  - Picks the correct topographic input data tiles into INPUT_ORO
  - Runs a batched job for Daint or Levante
  - Wrapper needs same Python environment (PYTHONPATH) as Extpar
  - run `WrapExtpar.py -h` for more infos
* Testsuite
  - Add new test for WrapExtpar.py on Daint and Levante
* Python
  - Use `raise` instead of `sys.exit(1)` to abort scripts


## 5.10.1

This is a minor release to fix cdo version on Levante
* CDO-version fixed to 2.0.5 because later versions exit with non-zero exit status for `cdo -V`

## 5.10

This is a minor release with fixes for Levante at DKRZ, a cleanup of MCH runscripts and a bugfix for in consistency_check
* CDO version change from 1.9.10 to 2.0.5 on Levante at DKRZ
   - Fields from Python-CDO scripts changed up to 10e-4!
   - Many (even large) adjustements of the tolerances in the testsuite
   - Result on other machines remain unchanged
* Fix for array out-of-bounds with special points outside of target domain
* Abort Extpar in consistency check if compiler not GCC
* Cleanup and update runscripts for models run at MeteoSwiss

## 5.9.1

This is a minor release to support Levante at DKRZ and drop the support for Intel
* Support GCC on Levante
* Remove configure-wrappers for Mistral
* Remove Intel references in testsuite

## 5.9

This is an intermediate release with OpenMP optimizations, a bugfix for topography and and some cleanup
* OpenMP optimizations for domains crossing date-line in the following parts
   - extpar_landuse_to_buffer
   - extpar_topo_to_buffer
* Update runscripts
   - ETH domain with Merit topography
   - EU-Cordex domain with Globe topography
* Cleanup
   - Remove all leftovers from old build-system
   - Remove output for vertices
      - Around 10% speedup

* Bugfix: Take last row at dateline into account
   - Last raw data row in mo_agg_topo was ignored
   - Changes the following fields (only in testcase mpim/icon_r2b4)
      - SSO_GAMMA
      - SSO_OROMAX
      - SSO_OROMIN
      - SSO_SIGMA
      - SSO_STDH
      - SSO_THETA
      - T_CL
      - Z0

## 5.8

This is an intermediate release with changes for the upgrade of Piz Daint, support for Merit topography for COSMO and a revised algorithm for SGSL processing
* Daint upgrade
   - export PMI_NO_PREINITIALIZE=1 to avoid unwanted prints from CDO
   - New paths for Python virtual environment:
      - Daint: /project/g110/extpar/venv_daint
      - Tsa: /project/g110/extpar/venv_tsa
* Merit for COSMO
   - Enable itopo_type=3 for COSMO grid
   - Testing performed by Christian Steger from Hymet-group at ETH-IAC
* Revised algorithm for SGSL preprocessing
   - Infer values of boundary points for SGSL
   - Deactivate SGSL for ASTER due to inconsistencies. For detailed information see the [PR :material-open-in-new:](https://github.com/C2SM-RCM/extpar/pull/278){:target="_blank"}.

## 5.7.4

This is a minor release with an adaption in the SSO-computation and replacement of ksh with bash
* Revised SSO for ICON
   - Experiments at DWD showed an improvement of model results with different SSO-thresholds, see this [wiki entry :material-open-in-new:](https://github.com/C2SM-RCM/extpar/wiki/Reducing-SSO-STDH-limit){:target="_blank"} for more information.
   - Lowering of threshold from 10.0 to 1.0 changes the following fields
      - SSO_GAMMA
      - SSO_OROMAX
      - SSO_OROMIN
      - SSO_SIGMA
      - SSO_STDH
      - SSO_THETA
      - topography_c (only in testcase mpim/icon_r2b4)
     
* Replace ksh with bash in runscripts to prepare for future linux-distributions

## 5.7.3

This is a minor release with two technical improvements.
* Python-CDO
   - Automatic detection if CDO contains thread-safe HDF5 library
* extpar_topo_to_buffer
   - Reduce memory usage for option *lsubstract_mean_slope=.TRUE.*

## 5.7.2

This is a minor release with an update for the Python environment on Mistral and refactored docs.
* Python on Mistral
  - Replace anaconda3/bleeding_edge with python3/unstable
* Docs
  - Remove outdated documentation
  - Make docs about testsuite and compilation more comprehensive
  - Quickstart for all supported machines

## 5.7.1

This is a minor release with two bugfixes, one for the build-system and one for extpar_topo_to_buffer.
* Build-system:
  - Install Python-CDO script in bin-folder directly in Makefile
  - Out-of-source build include Python-CDO scripts as well
 * extpar_topo_to_buffer:
   - Abort Extpar for itopo_type = 3

## 5.7

This is an intermediate release that introduces two Python-CDO scripts, modifications for reduced memory usage
for non-global grids and a bugfix for the CAMS-aersosol dataset.
* extpar_isa_to_buffer
   - Replace Fortran code with Python-CDO
   - Change in results for fields
        - Impervous Surface Area (ISA)
   - Read the [users guide](doc/user_and_implementation_manual.pdf) for detailed information about details of the implementation of *extpar_isa_to_buffer*
* exptar_ahf_to_buffer
   - Replace Fortran code with Python-CDO
   - Change in results for fields
      - Antropogenic Heat Flux (AHF)
   - Read the [users guide](doc/user_and_implementation_manual.pdf) for detailed information about details of the implementation of *extpar_ahf_to_buffer*
* Reduce memory usage for Python-CDO
   - Automatic determination of the extent of the target grid
   - Use CDO operator **-sellonlat** to read subset of input data
* Bugfix for CAMS dataset
   - Remove assertion for GCC-compiler
   - Fix wrong dimensions in meta-data for the ICON-grid

## 5.6

This is an intermediate release that introduces a new topography dataset and the CAMS-aerosol climatologies, 
OpenMP support for CSCS-machines, enhanced testing on CSCS-machines, a script to extract the input-data from namelist, 
progress bars for logfiles of extpar_topo_to_buffer and consistent names for all logfiles.

* Merit/Rema topography
  - Set switch itopo_type = 3 to process Merit-Rema data
* CAMS aerosol climatology for ICON
  - Set switch iaot_type = 5 to process CAMS-aerosol data
  - Due to an unresolved bug only available for Intel compiler
* Testsuite
  - Remove COSMO-D2 test from testsuite
  - Enable landuse tests for all supported machines
  - All tests run on 12 OpenMP threads on Piz Daint and Tsa for COSMO and ICON
  - Introduce the script [extract_inputfiles_from_namelist.py](test/testsuite/bin/extract_inputfiles_from_namelist.py) for faster data access at CSCS. It is recommended to use this skript for your own Extpar runs as well.
* Logging
  - Progress bar (0% to 100%) for topography processing
  - Change logfile-name for extpar_consistency_check to extpar_consistency_check.log

## 5.5.1

This is a minor release that fixes the inconsistent usage of netCDF versions across Extpar and small documentation changes.

* Bugfix netCDF versions
   - Use netCDF version passed via *NETCDF_OUTPUT_FILETYPE* also for ICON grids (CDI-interface)
* Documentation changes
   - Put Jonas Jucker as source code administrator

## 5.5

This is an intermediate release that brings enhanced namelist parsing for the Python-CDO scripts, a new Python-CDO script *extpar_era_to_buffer.py* to replace
the former way of remapping ERA-climatologies using Icontools, a more sophisticated tolerance checker to allow specific roundoff for each test and variables, support NetCDF5 and new default NetCDF 4, fixes for high-resolution grid exceeding integer value range and some minor bugfixes for Piz Daint related to HDF5.  
* exptar_era_to_buffer  
   - 4 fields processed
      - Sea surface temperature (T_SEA)
      - 2m land temperature (T_2M_CLIM)
      - Geometric height (TOPO_CLIM)
      - Snow water equivalent (W_SNOW)
   - New namelist-parameter iera_type defines type of ERA input data used, either ERA-I or ERA-5
   - *extpar_consistency_check* checks for namelist *INPUT_ERA* to determine if ERA-climatologies come from Python-CDO or Icontools
   - Using *extpar_era_to_buffer.py* changes fields, a detailed [review of changes :material-open-in-new:](https://github.com/C2SM-RCM/extpar/wiki/Comparison-of-ERA-Interim-processing-using-ICON-REMAP-and-Python-CDO){:target="_blank"} was performed by Jürgen Helmert from DWD
      - W_SNOW
      - TOPO_CLIM
      - T_SEA
      - T_2M_CLIM
   - Read the [users guide](doc/user_and_implementation_manual.pdf) for detailed information about how *extpar_era_to_buffer* is integrated into the existing workflow
* Enhanced namelist parsing for Python-CDO
   - Line starting with ! ignored as expected from Fortran code
* Bugfixes for Piz Daint
   - -L option for all CDO commands
   - Disbable HDF5 file locking due to problems reading some input data
* Improved tolerance testing in testsuite
   - Tolerances can now be defined separate for each test and variable for example in [tolerance file](test/testsuite/data/dwd/icon_d2/tolerances)
* Support for NetCDF 5
  - NetCDF 4 replaces netCDF 3 as default
  - Value of environment variable *NETCDF_OUTPUT_FILETYPE* sets version: NETCDF3, NETCDF4 or NETCDF5
* Modified netCDF-interface functions to allow write of fields with dimesions exceeding default integer value range 
      
## 5.4.1

This is an intermediate release that brings two lradtopo-parameters for Icon, better user feedback for the shell-commands launched in the Python-scripts, a bugfix in exptar_albedo_to_buffer.py, a configure script for O3 (ETHZ) and small technical improvements to the Code.
* HORIZON and SKYVIEW fields for the Icon grid
   - 4 new namelist-parameter
      - radius -> defines the considered horizontal distance for the HORIZON field
      - min_circ_cov -> defines the level of detail of the search-algorithm for performance reasons
      - max_missing -> defines upper treshold for the allowed missingness at the boundary of the domain
      - itype_scaling -> choose the type of scaling for SKYVIEW to account for anisotropic behaviour of IR-radiation
    - Read the [users guide](doc/user_and_implementation_manual.pdf) for detailed information about the difference between the COSMO and the ICON implementation
    
* Refactor function *launch_shell* by using subprocess.PIPE, providing output even when command crashes
* Correct bug for ialb_type=1 or 2 during netcdf write
* Configure script for O3 at ETHZ, **not regularly tested with Jenkins**
* Change link to CDI-submodule, to allow access for people witout DKRZ account
* Split chained CDO-operators into two steps to prevent crashes on Piz  Daint

## 5.4

This is a major release that introduces a rewrite of 4 Extpar programmes in Python, a common git-LFS input data repository,
a new build-system, 2 additional landuse data sets, CDI-library for icon grids in consistency check, mmap-caching for consistency check for less memory usage, some small improvements in the Fortran code and some minor changes in the testsuite.

* Rewrite of 4 Extpar programmes in Python
   - Modules extpar_alb_to_buffer.py, extpar_cru_to_buffer.py, extpar_emiss_to_buffer.py and extpar_ndvi_to_buffer.py
   - Small changes of the fields compared to the former Fortran implementation due to different interpolation methods, especially at the coastlines
   - Fields changed:
      - NDVI, NDVI_MAX, NDVI_MRAT
      - ALB_SAT, ALB_DRY
      - ALB, ALUVD, ALNID
      - T_CL
      - EMISS_RAD
   - A [review :material-open-in-new:](https://github.com/C2SM-RCM/extpar/wiki/Review-of-fields-for-Extpar-Version-5.4){:target="_blank"} involving users from DWD, MCH, MPIM and ETH took place to ensure the correctness of all fields changed
   - All Python programmes read from the same namelist file *namelist.py* containing Python dictionaries for each Extpar program.
   - Support of the old and coarse data (it_cl_type = 2) in extpar_cru_to_buffer expires and is replaced the following:
      - it_cl_type = 2 aggregates the coarse data over sea **and** the fine data over land
      - it_cl_type = 1 aggregates the fine data over land, sea points are not considered
      - For aggregation of the coarse data over land and sea only, use Extpar 5.3 or older
      
   - Read the [users guide](doc/user_and_implementation_manual.pdf) for detailed information about the rewritten programmes.
   
* git-LFS input data repository
   - All input data that can be processed with Extpar is stored in a unified data repository [extpar-input-data :material-open-in-new:](https://gitlab.dkrz.de/extpar-data/extpar-input-data){:target="_blank"}
   - Move all useful scripts and informations from raw_data_tools to the data repository hosted at DKRZ.
   - Remove folder raw_data_tools from Extpar repository
   - Some fields are renamed for better understanding, so please check your runscripts to adapt the new names.
   - Location on CSCS: */store/c2sm/extpar_raw_data/linked_data*
   - Location on Mistral: */work/pd1167/extpar-input-data/linked_data*
* New build-system  
   - New system follows the configure/make/make install paradigm
   - Out-of-source build supported
   - 4 basic steps to compile Extpar into binaries:
      - Run configure.*your_machine*.*your_compiler*
      - source modules.env
      - make or make -j 4 (for faster compilation)
    - Kesch at CSCS is no longer supported
   
* Corine landuse data
   - Additional landuse data set covering Europe
   - Can only be used in combination with GLOBCOVER (i_landuse_data=1)
   - Set switch l_use_corine=.true. in namelist *lu_raw_data* to aggregate the new data set
   - The corine landuse data set is only tested on Mistral at DKRZ

* ECCI landuse data
   - Global landuse data set split in 6 tiles
   - Set switch i_landuse_data=5, ntiles_globcover=6 and ilookup_table_lu = 1 in namelist *lu_raw_data* to aggregate the new data set
   - The ECCI landuse data is only tested on Mistral at DKRZ
   
* Enhanced testsuite
   - Icon test for DWD for all compilers
   - Jenkins on Mistral, Tsa and Daint
   - Convert testsuite src-code from Python2 to Python3
   - Pep8-Coding style test for Python code
   - Allow round-off for certain fields in output
   - Copy all required files from namelistdir (icon grids, clim-fields and Python-files) through testsuite itself

* CDI library for icon grids
   - [CDI :material-open-in-new:](https://code.mpimet.mpg.de/projects/cdi){:target="_blank"} write routine replaces write_netcdf_icon_grid routine
   - Output of icon grids **always** involves CDI, output without CDI no longer supported
   - CDI contained as a git submodule inside the Extpar repository
   - See [compile_run](doc/README.compile_run.md) for instructions to clone Extpar from GitHub correctly
   
* Mmap-caching
   - allows run of Extpar on machines with only little memory
   - new logical parameter *l_use_array_cache = .true. * in namelist file *INPUT_CHECK* activates mmap-caching
   - Bitwise-identical with and without mmap-caching
   - Only supported and tested for GCC compiler

* Fortran Code changes
   - Remove all *filename_max* from INTENT(IN)
   - Output of *COSMO/ICON* netCDF-files in the buffer modules no longer supported
   - Remove all unused modules/programmes replaced by Python modules
   - Remaining code still needed in Fortan now contained in modules *mo_python_data.f90*,  
     *mo_python_routines.f90*, and *mo_python_tg_fields.f90*
   

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

### Build Mechanism

 * Added Options file for compiling on LCE with Intel compiler
 * Added Options files for compiling on Mistral with GCC, NAG, and Intel compilers
 * Added Options files for compiling on Mac OS with the GCC compiler
 * Updated Options file for compiling on o3 at ETH with PGI compiler
 * Update bin/gen_info.sh for compatibility with git

### Run scripts

 * Minor adaptations to run scripts to be compatible with new version
 * Addition of MPI ICON run script

### Testing

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
