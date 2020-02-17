#!/bin/ksh

run_sequential()
{
    script=$1
    set +e
    echo ">> Run ${script} ..."   >>  ${logfile}
    start=$(date +%s.%N)
    ./${script}
    rc=$?
    printf "   Return code: %i\n" $rc >> ${logfile}
    end=$(date +%s.%N)
    (( runtime = end - start ))
    if [[ $rc -eq "0" ]]; then

        echo "   SUCCESS ${1%% *}" >> ${logfile}
        echo "   execution time: $runtime s" >> ${logfile}
    else
        case $rc in
            127)
                echo "   ERROR ${1%% *}: command not found" >> ${logfile}
                ;;
            130)
                echo "   ERROR ${1%% *}: script terminated by Ctrl-C" >> ${logfile}
                ;;             
            *)
                echo "   ERROR ${1%% *}: fatal error - return code $rc" >> ${logfile}
                echo "   See logfile of ${script} for more detailed information" >> ${logfile}
                ;;
        esac

            echo "*****************************************"
            echo ""
            echo "Some Extpar executables did not terminate correctly!"
            echo "See ${logfile} for more information"
            echo ""
            echo "*****************************************"
            exit 
    fi

    set -e
}

# function to launch Extpar executables in parallel
run_parallel()
{
    script=$1
    (echo ">> Run ${script} ..."  ;  (time ./${script})> ./time_$script 2>&1 ; echo $? > ./exit_status_$script )&
}

# function to check exit status and print timings of parallel launched Extpar executables 
check_exit_status()
{
    script=$1
    typeset -n count=$2

    exit_status=`cat exit_status_$script`
    if [[ $exit_status -eq "0" ]]; then
        echo ""  >>  ${logfile}
        echo ">>${script} succesfully terminated in `cat time_$script`" >> ${logfile}
    else
        (( count++ ))
        echo "" >> ${logfile}
        case $exit_status in

            127)
                echo "   ERROR ${1%% *}: command not found" >> ${logfile}
                ;;
            130)
                echo "   ERROR ${1%% *}: script terminated by Ctrl-C" >> ${logfile}
                ;;             
            *)
                echo "   ERROR ${1%% *}: fatal error - return code $exit_status" >> ${logfile}
                echo "   See logfile of ${script} for more detailed information" >> ${logfile}
                ;;
        esac
    fi
}
