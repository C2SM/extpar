#!/bin/bash
#
# Script for creating info_extpar.f90 from a template replacing the old
# way to be more compatible with Git.
#
# Author  : Luis Kornblueh, Max Planck Institute for Meteorology
#           luis.kornblueh@mpimet.mpg.de 
#
# History :
#   2018-03-23  First release
#   2018-03-23  Updated to work properly with new build system and out-of-source build
#_____________________________________________________________________________________________
#
set -eu

# parse command line arguments
fconfig="$1"
srcdir="$2"
builddir="${3-${srcdir}}"

extpar_dir=${2%/src*}

#_____________________________________________________________________________________________
#
# set some globals (need to be manually updated)
code_name="extpar"
code_version=$(cd ${extpar_dir} && git describe --tags | awk -F'-' '{print $1}')

#_____________________________________________________________________________________________
#
# check if file is available
if [[ ! -f "${srcdir}/info_extpar.f90.in" ]]
then
  echo "ERROR: could not find ${srcdir}/info_extpar.f90.in file" 1>&2
  exit 1
fi

#_____________________________________________________________________________________________
#
# get host and date
creation_host=$(hostname | sed 's/[0-9]*$//g')
creation_date=$(date -u +"%F %T")
creation_user=$(whoami)

#_____________________________________________________________________________________________
#
# get compiler information
compiler=$(grep 'Compiler command' $fconfig | cut -d":" -f2)
cversion=$(grep 'Compiler version' $fconfig | cut -d":" -f2)
cincludes=$(grep 'Compiler includes' $fconfig | cut -d":" -f2)
cflags=$(grep 'Compiler flags' $fconfig | cut -d":" -f2)
linker=$(grep 'Linker command' $fconfig | cut -d":" -f2)
lversion=$(grep 'Linker version' $fconfig | cut -d":" -f2)
lflags=$(grep 'Linker flags' $fconfig | cut -d":" -f2)
llibraries=$(grep 'Linker libraries' $fconfig | cut -d":" -f2 |cut -c1-200)
string=$(grep 'Compiler *:' $fconfig | cut -d":" -f2 | cut -d" " -f3-)

#_____________________________________________________________________________________________
#
# information related to version control system
inside_git_repository=""
cd ${extpar_dir}
inside_git_repository=$(git rev-parse --is-inside-work-tree 2> /dev/null)

if [[ "$inside_git_repository" = "true" ]]
then
    # this is a git clone copy
    INFO_PackageName="${code_name}-${code_version}"
    INFO_RepositoryURL=$(git --git-dir ./.git remote get-url origin)
    INFO_LastCommitDate=$(git log -1 --format=%cd --date=iso | awk -F" " '{print $1, $2}')
    INFO_RevisionHash=$(git rev-parse HEAD)
    INFO_CodeIsModified=$([[ $(git status 2> /dev/null | tail -n1) != "nothing to commit, working tree clean" ]] && echo "modified" || echo "clean")
else
    # this directory is not under version control
    INFO_PackageName="${code_name}-${code_version}"
    INFO_RepositoryURL="(missing)"
    INFO_LastCommitDate="(missing)"
    INFO_RevisionHash="(missing)"
    INFO_CodeIsModified="unknown"
fi

cd - > /dev/null 2>&1

#_____________________________________________________________________________________________
#
# information related to compilation
INFO_CompilerVersion="${cversion/ gnu/GCC}"
INFO_CompiledBy="${creation_user}"
INFO_CompileTime="${creation_date}"
INFO_CompileMachine="${creation_host}"

#_____________________________________________________________________________________________
#
# information which will have to be defined at runtime
INFO_StartTime=""
INFO_BinaryName=""


#_____________________________________________________________________________________________
#
# replace the placeholders in info_extpar.f90 with the correct strings
sed -e "s|\(INFO_PackageName *= *\)'.*'|\1'""${INFO_PackageName}""'|g
        s|\(INFO_RepositoryURL *= *\)'.*'|\1'""${INFO_RepositoryURL}""'|g
        s|\(INFO_LastCommitDate *= *\)'.*'|\1'""${INFO_LastCommitDate}""'|g
        s|\(INFO_RevisionHash *= *\)'.*'|\1'""${INFO_RevisionHash}""'|g
        s|\(INFO_CodeIsModified *= *\)'.*'|\1'""${INFO_CodeIsModified}""'|g
        s|\(INFO_CompilerVersion *= *\)'.*'|\1'""${INFO_CompilerVersion}""'|g
        s|\(INFO_CompiledBy *= *\)'.*'|\1'""${INFO_CompiledBy}""'|g
        s|\(INFO_CompileTime *= *\)'.*'|\1'""${INFO_CompileTime}""'|g
        s|\(INFO_CompileMachine *= *\)'.*'|\1'""${INFO_CompileMachine}""'|g
        s|\(INFO_StartTime *= *\)'.*'|\1'""${INFO_StartTime}""'|g
        s|\(INFO_BinaryName *= *\)'.*'|\1'""${INFO_BinaryName}""'|g" ${srcdir}/info_extpar.f90.in > ${builddir}/info_extpar.f90

exit 0



