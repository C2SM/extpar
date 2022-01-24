#!/bin/bash

# This script runs the Extpar testsuite 
#

#######################################################################################
# Function section
# Function to launch slurm jobs and track when job is finished
function launch_job {
      local script=$1
      local timeout=$2

              # check sanity of arguments
      test -f "${script}" || exitError 7201 ${LINENO} "cannot find script ${script}"
      if [ -n "${timeout}" ] ; then
        echo "${timeout}" | grep '^[0-9][0-9]*$' 2>&1 > /dev/null
        if [ $? -ne 0 ] ; then
          exitError 7203 ${LINENO} "timeout is not a number"
        fi
      fi
      # get out/err of SLURM job
      local out=`grep '^\#SBATCH --output=' ${script} | sed 's/.*output=//g'`

      # submit SLURM job
      local res=`sbatch ${script}`
      if [ $? -ne 0 ] ; then
        exitError 7205 ${LINENO} "problem submitting SLURM batch job"
      fi
      echo "${res}" | grep "^Submitted batch job [0-9][0-9]*$" || exitError 7206 ${LINENO} "problem determining job ID of SLURM job"
      local jobid=`echo "${res}" | sed  's/^Submitted batch job //g'`
      test -n "${jobid}" || exitError 7207 ${LINENO} "problem determining job ID of SLURM job"
      # wait until job has finished (or maximum sleep time has been reached)
      if [ -n "${timeout}" ] ; then
        local secs=0
        local inc=2
        local job_status="UNKNOWN"
        while [ $secs -lt $timeout ] ; do
             echo "...waiting ${inc}s for SLURM job ${jobid} to finish (status=${job_status})"
             sleep ${inc}
             secs=$[$secs+${inc}]
             inc=60
             squeue_out=`squeue -o "%.20i %.20u %T" -h -j "${jobid}" 2>/dev/null`
             echo "${squeue_out}" | grep "^ *${jobid} " &> /dev/null
             if [ $? -eq 1 ] ; then
               break
             fi
             job_status=`echo ${squeue_out} | sed 's/.* //g'`
        done
      fi
      # make sure that job has finished
      squeue_out=`squeue -o "%.20i %.20u %T" -h -j "${jobid}" 2>/dev/null`
      echo "${squeue_out}" | grep "^ *${jobid} " &> /dev/null
      if [ $? -eq 0 ] ; then
          exitError 7207 ${LINENO} "batch job ${script} with ID ${jobid} on host ${slave} did not finish"
      fi
}
# Error function
exitError()
{
        echo "ERROR $1: $2" 1>&2
        echo "ERROR     LOCATION=$0" 1>&2
        exit $1
}

# Code body

case "$(hostname)" in
    daint*)
        host=daint
        module load cray-python
        ;;
    tsa*)
        host=tsa
        module load python/3.7.4
        ;;
    mlogin*)
        host=mistral
        source /etc/profile.d/mistral.sh
        module unload cdo
        module load cdo
	;;
esac

cd test/testsuite

# Get the input data
cd data
./get_data.sh
cd ..

# Extract data files from namelists
./bin/extract_inputfiles_from_namelist.py


# Copy the executables
cp ../../bin/* bin

if [[ "$host" == "daint" || "$host" == "tsa" ]]; then
    echo "Running transfer script"
    script="./submit.${host}.transfer.sh"
    test -f ${script} || exitError 1260 "submit script ${script} does not exist" 
    launch_job ${script} 7200
    if [ $? -ne 0 ] ; then
      exitError 1251 ${LINENO} "problem launching SLURM job ${script}"
      cat transfer.log
    fi
    echo "Finished with transfer script"
    cat transfer.log
fi

if [ "$compiler" = "intel" ]; then
  script="./submit.mistral.intel.sh"
else
  script="./submit.${host}.sh"
fi
test -f ${script} || exitError 1260 "submit script ${script} does not exist" 

echo "Running submit script"
launch_job ${script} 7200
if [ $? -ne 0 ] ; then
  exitError 1251 ${LINENO} "problem launching SLURM job ${script}"
fi
echo "Finished with submit script"

# echo output to stdout
test -f testsuite.out || exitError 1261 "output file testsuite.out does not exist"
echo "=== testsuite.out BEGIN ==="
cat testsuite.out | /bin/sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g"
echo "=== testsuite.out END ==="

# check result of testsuite
grep RESULT testsuite.out | egrep 'FAIL|CRASH' > /dev/null
if [ $? -eq 0 ] ; then
  exitError 1271 "testsuite did not complete successfully"
fi
