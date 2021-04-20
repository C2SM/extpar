AC_PREREQ([2.69])
AC_INIT([extpar], [5.3.0])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/extpar_aot_to_buffer.f90])

AC_PROG_CC
AC_LANG([C])
dnl Supporting ISO C89 means that the following header files are available:
dnl stdlib.h, stdio.h, math.h, time.h, ctype.h, string.h, signal.h, limits.h
dnl errno.h, ctype.h assert.h.
AS_VAR_IF([ac_cv_prog_cc_c89], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C89])])

dnl Supporting ISO C99 means that the following header files are available:
dnl stdbool.h.
AC_PROG_CC_C99
AS_VAR_IF([ac_cv_prog_cc_c99], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C99])])

ACX_COMPILER_CC_VENDOR_SIMPLE
ACX_COMPILER_CC_VERSION_SIMPLE

AC_PROG_FC
AC_LANG([Fortran])
ACX_FC_PP_SRCEXT([f90])
ACX_COMPILER_FC_VENDOR_SIMPLE
AC_SUBST([FC_VENDOR], [$acx_cv_fc_compiler_vendor])
ACX_COMPILER_FC_VERSION_SIMPLE
AC_SUBST([FC_VERSION], [$acx_cv_fc_compiler_version])

dnl ACX_FC_LINE_LENGTH([unlimited])

ACX_FC_MODULE_NAMING(
  [AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])
   AC_SUBST([FCMODUC], ["$acx_cv_fc_module_naming_upper"])])

ACX_FC_MODULE_IN_FLAG([AC_SUBST([FCMODINC], ["$acx_cv_fc_module_in_flag"])])
ACX_FC_MODULE_OUT_FLAG([AC_SUBST([FCMODOUT], ["$acx_cv_fc_module_out_flag"])])

ACX_FC_INCLUDE_FLAG([AC_SUBST([FCINCFLAG], ["$acx_cv_fc_ftn_include_flag"])])
ACX_FC_INCLUDE_ORDER(
  [AC_SUBST([FCINCORDER], ["$acx_cv_fc_ftn_include_order"])],
  [AC_SUBST([FCINCORDER], ['src,flg'])])

ACX_FC_INCLUDE_FLAG_PP(
  [AC_SUBST([FCINCFLAG_PP], ["$acx_cv_fc_pp_include_flag"])])
ACX_FC_INCLUDE_ORDER_PP(
  [AC_SUBST([FCINCORDER_PP], ["$acx_cv_fc_pp_include_order"])],
  [AC_SUBST([FCINCORDER_PP], ['inc,flg'])])

ACX_LANG_MACRO_FLAG([AC_SUBST([FCDEF_PP], ["$acx_cv_fc_macro_flag"])])

dnl Additional macro definitions to aid the dependency generator,
dnl e.g. compiler intrinsic macros:
AC_SUBST([DEPGEN_FCFLAGS], [''])

dnl A list of flags enabling location of the include files of the enabled
dnl bundled libraries:
AC_SUBST([BUNDLED_FCFLAGS], [''])

dnl A list of library files (i.e. *.a files) of the enabled bundled libraries:
AC_SUBST([BUNDLED_LIBFILES], [''])

dnl Optionally enable silent rules
AC_ARG_ENABLE([silent-rules],
  [AC_HELP_STRING([--enable-silent-rules],
     [less verbose build output (undo: "make V=1") @<:@default=no@:>@])], [],
  [enable_silent_rules=no])

AS_VAR_IF([enable_silent_rules], [yes],
  [AC_SUBST([DEFAULT_VERBOSITY], [0])],
  [AC_SUBST([DEFAULT_VERBOSITY], [1])])

dnl Generate RPATH flags if needed:
AC_ARG_ENABLE([rpaths],
  [AC_HELP_STRING([--enable-rpaths],
     [add directories specified with -L flags in LDFLAGS and LIBS to the
runtime library search paths (RPATH) @<:@default=yes@:>@])], [],
  [enable_rpaths=yes])

AC_SUBST([RPATH_LDFLAGS], [''])
AS_VAR_IF([enable_rpaths], [yes],
  [ACX_SHLIB_FC_RPATH_FLAG
   ASX_EXTRACT_ARGS([acx_L_args], ["$LDFLAGS $LIBS"], ['-L@<:@ @:>@*'])
   for acx_L_arg in $acx_L_args; do
     ASX_VAR_APPEND_UNIQ([RPATH_LDFLAGS],
       ["$acx_cv_fc_rpath_flag$acx_L_arg"], [' '])
   done
   ACX_SHLIB_RPATH_FLAGS_CHECK([$RPATH_LDFLAGS])])

dnl Check for OpenMP
AC_ARG_ENABLE([openmp],
  [AC_HELP_STRING([--enable-openmp],
     [enable OpenMP support @<:@default=no@:>@])], [],
  [enable_openmp=no])

AS_VAR_IF([enable_openmp], [yes],
  [ACX_LANG_OPENMP_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openmp_flag"])
      ACX_LANG_MACRO_CHECK_VALUE([_OPENMP],
         [201511 201307 201107 200805 200505 200011 199911 199710])
      AS_VAR_IF([acx_macro_value], [unknown],
        [AC_MSG_WARN([unable to detect supported OpenMP version])
         AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FCDEF_PP}_OPENMP"])],
        [AS_VAR_APPEND([DEPGEN_FCFLAGS],
           [" ${FCDEF_PP}_OPENMP=$acx_macro_value"])])])],
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENMP])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenMP support by default: set dnl
FCFLAGS accordingly to disable it])])])

dnl Check for NetCDF Fortran
AC_ARG_WITH([netcdf-fortran],
  [AC_HELP_STRING([--with-netcdf-fortran=<yes|directory>],
     [location of NetCDF Fortran library (lib and include subdirs)])],
  [AS_VAR_IF([withval], [no],
     [AC_MSG_ERROR([cannot build the package without NetCDF Fortran dnl
library])])],
  [with_netcdf_fortran=yes])
AC_ARG_VAR([NETCDF_FCFLAGS],
  [exact Fortran compiler flags enabling NetCDF Fortran])
AC_ARG_VAR([NETCDF_FCLIBS],
  [exact linker flags enabling NetCDF Fortran when linking with Fortran ]dnl
[compiler])

acx_inc_search_args=
acx_lib_search_args=
AS_IF([test x"$with_netcdf_fortran" != xyes],
  [acx_inc_search_args="${FCMODINC}$with_netcdf_fortran/include"
   acx_lib_search_args="-L$with_netcdf_fortran/lib"
   AS_VAR_IF([enable_rpaths], [yes],
     [AS_VAR_APPEND([acx_lib_search_args],
        [" ${acx_cv_fc_rpath_flag}$with_netcdf_fortran/lib"])])])

acx_save_FCFLAGS=$FCFLAGS
NETCDF_FCFLAGS=${NETCDF_FCFLAGS-$acx_inc_search_args}
FCFLAGS="$NETCDF_FCFLAGS $acx_save_FCFLAGS"
ACX_FC_MODULE_CHECK([netcdf],
  [ACX_LANG_LIB_SEARCH([NETCDF_FCLIBS], [nf90_open],
     [ASX_PREPEND_LDFLAGS([$acx_lib_search_args],
        [-lnetcdff], [-lnetcdff -lnetcdf])],
     [NETCDF_FCLIBS=$acx_cv_fc_lib_func_NF90_OPEN
      acx_have_netcdf=yes], [],
     [AC_LANG_PROGRAM([],
[[      use netcdf, only: nf90_open, NF90_NOWRITE
      implicit none
      integer :: status, ncid
      status = nf90_open('conftest.nc', NF90_NOWRITE, ncid)]])])])
FCFLAGS=$acx_save_FCFLAGS

dnl Check for CDI
AC_ARG_WITH([cdi],
  [AC_HELP_STRING([--with-cdi=<yes|directory|bundled>],
     [location of CDI library (lib and include subdirs)])],
  [AS_VAR_IF([withval], [no],
     [AC_MSG_ERROR([cannot build the package without CDI library])])],
  [with_cdi=bundled])
AC_ARG_VAR([CDI_FCFLAGS],
  [exact Fortran compiler flags enabling CDI])
AC_ARG_VAR([CDI_FCLIBS],
  [exact linker flags enabling CDI when linking with Fortran compiler])

AS_VAR_IF([with_cdi], [bundled],
  [AS_IF([test ! -r "$srcdir/bundled/cdi/src/cdi.h"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of CDI dnl
library (src/cdi.h) in '$srcdir/bundled/cdi': provide an external version dnl
of CDI (--with-cdi)])])
   AC_DISABLE_OPTION_CHECKING
   ACX_CONFIG_SUBDIR([bundled/cdi],
     [[F77=*|CXX=*|LDFLAGS=*|LIBS=*],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([static])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([shared])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([maintainer-mode])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([silent-rules])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([iso-c-interface])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([cf-interface])],
      [ACX_CONFIG_SUBDIR_PATTERN_ENABLE([cdi-app])]],
     [[--enable-static], [--disable-shared], [--disable-maintainer-mode],
      [--enable-silent-rules="$enable_silent_rules"],
      [--enable-iso-c-interface], [--disable-cf-interface],
      [--disable-cdi-app],
dnl Disable unwanted extra dependencies.
      [acx_cv_have_nc4hdf5=no],
dnl Disable the compilers we do not need.
      [F77=no],
dnl Note that we need to set CXX to 'no' not only to disable the redundant
dnl checks but also to avoid the inconsistency in the generation of the
dnl 'postdeps' variable for Fortran (!!!) compiler inside the *.la files. We
dnl care about those files because they are used when linking the test suite of
dnl CDI (i.e. when running 'make check'). If you want to know more, take a look
dnl at the macro _LT_LANG_FC_CONFIG in libtool.m4 and note that it does not
dnl initialize the variable 'output_verbose_link_cmd' as it should but uses the
dnl value that is set for the C++ compiler. Therefore, by disabling CXX we get
dnl the same behaviour (still questionable but satisfactory for us) as if CDI
dnl didn't use the C++ compiler at all.
      [CXX=no],
dnl We might have already updated LDFLAGS and LIBS and we want CDI to use them:
      [LDFLAGS="$LDFLAGS"], [LIBS="$LIBS"]],
dnl Output the help message of CDI's configure script when this configure script
dnl is called with '--help=recursive':
[yes])
   CDI_FCFLAGS=; CDI_FCLIBS=
   AS_VAR_APPEND([BUNDLED_FCFLAGS], ["${FCMODINC}bundled/cdi/src "])
   BUNDLED_LIBFILES="bundled/cdi/src/.libs/libcdi_f2003.a dnl
bundled/cdi/src/.libs/libcdi.a $BUNDLED_LIBFILES"],
  [acx_inc_search_args=
   acx_lib_search_args=
   AS_IF([test x"$with_cdi" != xyes && test x"$with_cdi" != xbundled],
     [acx_inc_search_args="${FCMODINC}$with_cdi/include"
      acx_lib_search_args="-L$with_cdi/lib"
      AS_VAR_IF([enable_rpaths], [yes],
        [AS_VAR_APPEND([acx_lib_search_args],
           [" ${acx_cv_fc_rpath_flag}$with_cdi/lib"])])])
   acx_save_FCFLAGS=$FCFLAGS
   CDI_FCFLAGS=${CDI_FCFLAGS-$acx_inc_search_args}
   FCFLAGS="$CDI_FCFLAGS $acx_save_FCFLAGS"
   ACX_FC_MODULE_CHECK([mo_cdi],
     [ACX_LANG_LIB_SEARCH([CDI_FCLIBS], [vlistDefVarTimetype],
        [ASX_PREPEND_LDFLAGS([$acx_lib_search_args],
           [-lcdi_f2003], [-lcdi_f2003 -lcdi])],
        [CDI_FCLIBS=$acx_cv_fc_lib_func_VLISTDEFVARTIMETYPE
         acx_have_netcdf=yes], [],
        [AC_LANG_PROGRAM([],
[[      use mo_cdi, only: vlistDefVarTimetype
      implicit none
      integer :: a, b, c
      call vlistDefVarTimetype(a, b, c)]])])])
   FCFLAGS=$acx_save_FCFLAGS
dnl Check for file formats supported by CDI.
   acx_save_FCFLAGS=$FCFLAGS; FCFLAGS="$CDI_FCFLAGS $acx_save_FCFLAGS"
   acx_save_LIBS=$LIBS; LIBS="$CDI_FCLIBS $acx_save_LIBS"
   acx_cdi_formats="grb:CDI_FILETYPE_GRB grb2:CDI_FILETYPE_GRB2
                    nc:CDI_FILETYPE_NC nc2:CDI_FILETYPE_NC2
                    nc4:CDI_FILETYPE_NC4 nc4c:CDI_FILETYPE_NC4C
                    nc5:CDI_FILETYPE_NC5 srv:CDI_FILETYPE_SRV
                    ext:CDI_FILETYPE_EXT ieg:CDI_FILETYPE_IEG"
   for acx_cdi_format in $acx_cdi_formats; do
     acx_cdi_format_name=`echo "$acx_cdi_format" | cut -d: -f1`
     acx_cdi_format_ID=`echo "$acx_cdi_format" | cut -d: -f2`
     acx_cache_var_name="acx_cv_cdi_have_$acx_cdi_format_name"
     AC_CACHE_CHECK([whether CDI supports file format '$acx_cdi_format_name'],
       [$acx_cache_var_name],
       [AS_VAR_SET([$acx_cache_var_name], [unknown])
        AS_VAR_IF([cross_compiling], [no],
          [AC_LINK_IFELSE([AC_LANG_PROGRAM([],
[[      use mo_cdi, only: cdiHaveFiletype, $acx_cdi_format_ID
      write(*, "(i0)") cdiHaveFiletype($acx_cdi_format_ID)]])],
             [acx_exec_result=`./conftest$ac_exeext 2>/dev/null`
              AS_IF([test $? -eq 0],
                [AS_IF([test $acx_exec_result -eq 1],
                   [AS_VAR_SET([$acx_cache_var_name], [yes])],
                   [AS_VAR_SET([$acx_cache_var_name], [no])])])])])])
   done
dnl Here we can check whether the format is supported and do something.
dnl For example, we can check whether NC is supported and emit a warning:
dnl    AS_IF([test x"$acx_cv_cdi_have_nc" != xyes],
dnl      [AC_MSG_WARN([external CDI does not support file format 'nc'])])
   FCFLAGS=$acx_save_FCFLAGS; LIBS=$acx_save_LIBS])

dnl The following is called at the very end, after the configure script of CDI
dnl is finished. Now we need to know what additional LDFLAGS and LIBS we need
dnl to be able to link to libcdi_f2003.a and libcdi.a.
AC_CONFIG_COMMANDS_PRE(
  [AS_VAR_IF([with_cdi], [bundled],
dnl First, we run config.status of CDI and retrieve LDFLAGS and LIBS collected
dnl by the configure script of CDI:
     [ACX_CONFIG_SUBDIR_VAR([acx_LDFLAGS], [bundled/cdi], [@LDFLAGS@])
      ACX_CONFIG_SUBDIR_VAR([acx_LIBS], [bundled/cdi], [@LIBS@])
dnl We passed all LDFLAGS and LIBS collected by this script to the configure
dnl script of CDI and we have not modified them since then. Therefore, it is
dnl safe to assume that acx_LDFLAGS and acx_LIBS have all the flags we need
dnl plus the flags collected by the configure script of CDI:
      LDFLAGS=$acx_LDFLAGS
      LIBS=$acx_LIBS
dnl We might need to add additional RPATH flags:
      AS_VAR_IF([enable_rpaths], [yes],
        [ASX_EXTRACT_ARGS([acx_L_args],
           ["$acx_LDFLAGS $acx_LIBS"], ['-L@<:@ @:>@*'])
         for acx_L_arg in $acx_L_args; do
           ASX_VAR_APPEND_UNIQ([RPATH_LDFLAGS],
             ["$acx_cv_fc_rpath_flag$acx_L_arg"], [' '])
         done])
dnl We might need to link to pthread library, which might require additional
dnl flags stored in PTHREAD_CFLAGS (we assume that the Fortran compiler
dnl understands the flags too):
      ACX_CONFIG_SUBDIR_VAR([acx_PTHREAD_CFLAGS],
        [bundled/cdi], [' @PTHREAD_CFLAGS@'])
dnl NAG compiler does not understand '-pthread' flag but it understands the flag
dnl '-Wl,-pthread', which does what we need here:
      AS_VAR_IF([FC_VENDOR], [nag],
        [acx_PTHREAD_CFLAGS=`AS_ECHO(["$acx_PTHREAD_CFLAGS"]) | dnl
sed 's/ -pthread/ -Wl,-pthread/g'`])
      AS_VAR_APPEND([LDFLAGS], ["$acx_PTHREAD_CFLAGS"])])])

dnl Check for PYTHON
AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
ACX_PROG_SEARCH([PYTHON], [python python3],
  [acx_python_try="$acx_candidate $srcdir/mkhelper/depgen.py -h dnl
>&AS_MESSAGE_LOG_FD"
   _AC_DO_VAR([acx_python_try])],
  [PYTHON=$acx_candidate])

AC_PROG_INSTALL
AC_PROG_MKDIR_P

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile])

dnl An additional step to generate dummy dependencies:
AC_CONFIG_COMMANDS([depfiles],
  [${MAKE-make} dummy-depend >/dev/null 2>&AS_MESSAGE_LOG_FD
   sleep 1
   touch Makefile])

AC_OUTPUT
