# Current Limitations

The EXTPAR software is subject to several limitations:

-   The ASTER domain can only be used from 60°N to 60°S. Be aware that an
    additional border of several gridpoints is needed if the
    topographically corrected parameters are desired. If the ASTER
    domain is exceeded a warning message is printed and the program
    *extpar_topo_to_buffer* is aborted.

-   The ASTER data shows some deficits, which are listed below:

    -   Beyond 60 degrees north and south, the ASTER raw data set
        features several areas where no value is available e.g., over
        Finland (private communication with HIRLAM).

    -   Some bogus regions may appear in complex topography. One of
        these regions is located near Grindelwald in Switzerland.

    -   The ASTER data are subject to artefacts of the satellite
        fly-over bands. Discontinuities can be spotted at the borders of
        such bands. In high latitudes these bands are better visible
        than in the low latitudes.

    -   As the correction of these deficits are time consuming no effort
        has been expended to remove these.

    -   The ASTER data might be subject to a shift of a half gridpoint
        (15 meters) in both directions.

-   There is no 3 km filtered ASTER or MERIT/REMA data set to derive the
    subgrid scale orography (`SSO`) parameters and the roughness length
    (`z0`) for ASTER.

-   The HWSD raw data is in a test phase. Furthermore a new version of
    Int2lm and TERRA is needed to make use of these data sets.

-   The special points are only tested for the COSMO grid. Also it is
    not possible to use these corrections if the soil raw data set is
    HWSD.

-   Array-caching in the consistency_check is only supported for GCC
    compiler.

-   CAMS aersosl data `iaot_type = 5` is only supported for Intel compiler.
