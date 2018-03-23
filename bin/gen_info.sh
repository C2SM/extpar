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
#
set -eu


# set some globals (need to be manually updated)
code_name="extpar"
code_version=$(git describe --tags $(git rev-list --tags --max-count=1) | cut -d"v" -f2)

# parse command line arguments
fconfig="$1"
srcdir="$2"

# check if file is available
if [ ! -f "${srcdir}/info_extpar.f90.in" ] ; then
  echo "ERROR: could not find ${srcdir}/info_extpar.f90.in file" 1>&2
  exit 1
fi

# get host and date
creation_host=$(hostname | sed 's/[0-9]*$//g')
creation_date=$(date -u +"%F %T")
creation_user=$(whoami)

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

# information related to version control system
if [[ ! -d .git ]]
then
    # this directory is not under version control
    INFO_PackageName="${code_name}-${code_version}"
    INFO_RepositoryURL="(missing)"
    INFO_LastCommitDate="(missing)"
    INFO_RevisionHash="(missing)"
    INFO_CodeIsModified="unknown"
else
    # this is a git clone copy
    INFO_PackageName="${code_name}-${code_version}"
    INFO_RepositoryURL=$(git --git-dir .git remote -v | grep fetch | awk -F" " 'NR==1{print $2}')
    INFO_CommitDate=$(git log -1 --format=%cd --date=iso | awk -F" " '{print $1, $2}')
    INFO_RevisionHash=$(git rev-parse HEAD)
    INFO_CodeIsModified=$([[ $(git status 2> /dev/null | tail -n1) != "nothing to commit, working tree clean" ]] && echo "modified" || echo "clean")
fi

# information related to compilation
INFO_CompilerCall="${compiler} ${cflags}"
INFO_CompilerVersion="${cversion}"
INFO_DefinedMacros="${cincludes}"
INFO_UndefinedMacros="(missing)"
INFO_DebugOptions="(missing)"
INFO_LinkOptions="${linker} ${llibraries}"
INFO_CompiledBy="${creation_user}"
INFO_CompileTime="${creation_date}"
INFO_CompileMachine="${creation_host}"

# information which will have to be defined at runtime
INFO_StartTime=""
INFO_BinaryName=""
INFO_RunMachine=""

# replace the placeholders in info_int2lm.f90 with the correct strings
sed -e "s|\(INFO_PackageName *= *\)'.*'|\1'""${INFO_PackageName}""'|g
        s|\(INFO_RepositoryURL *= *\)'.*'|\1'""${INFO_RpositoryURL}""'|g
        s|\(INFO_LastCommitDate *= *\)'.*'|\1'""${INFO_LastCommitDate}""'|g
        s|\(INFO_RevisionHash *= *\)'.*'|\1'""${INFO_RevisionHash}""'|g
        s|\(INFO_CodeIsModified *= *\)'.*'|\1'""${INFO_CodeIsModified}""'|g
        s|\(INFO_CompilerCall *= *\)'.*'|\1'""${INFO_CompilerCall}""'|g
        s|\(INFO_CompilerVersion *= *\)'.*'|\1'""${INFO_CompilerVersion}""'|g
        s|\(INFO_DefinedMacros *= *\)'.*'|\1'""${INFO_DefinedMacros}""'|g
        s|\(INFO_UndefinedMacros *= *\)'.*'|\1'""${INFO_UndefinedMacros}""'|g
        s|\(INFO_DebugOptions *= *\)'.*'|\1'""${INFO_DebugOptions}""'|g
        s|\(INFO_LinkOptions *= *\)'.*'|\1'""${INFO_LinkOptions}""'|g
        s|\(INFO_CompiledBy *= *\)'.*'|\1'""${INFO_CompiledBy}""'|g
        s|\(INFO_CompileTime *= *\)'.*'|\1'""${INFO_CompileTime}""'|g
        s|\(INFO_CompileMachine *= *\)'.*'|\1'""${INFO_CompileMachine}""'|g
        s|\(INFO_StartTime *= *\)'.*'|\1'""${INFO_StartTime}""'|g
        s|\(INFO_BinaryName *= *\)'.*'|\1'""${INFO_BinaryName}""'|g
        s|\(INFO_RunMachine *= *\)'.*'|\1'""${INFO_RunMachine}""'|g" ${srcdir}/info_extpar.f90.in > ${srcdir}/info_extpar.f90 

exit 0



