#!/bin/bash
case "$(hostname)" in
    *co2* | *iacdipl-7*)
        set -e

        MODE=${1:-'NORMAL'}
        if [[ ${MODE} == "DEBUG" ]]; then
            CURRENT_DATETIME=$(date --iso-8601=seconds)
            HASH=$(echo -n "${CURRENT_DATETIME}" | sha256sum | awk '{print $1}')
            mkdir -p /net/co2/c2sm-services/extpar/${HASH}

            podman run -e OMP_NUM_THREADS=16 -v /net/co2/c2sm-data/extpar-input-data:/data -v /net/co2/c2sm-services/extpar/${HASH}:/artifacts extpar:$ghprbPullId bash -c "/workspace/test/jenkins/test_docker.sh ${MODE}" || (podman image rm -f extpar:$ghprbPullId && exit 1)
        else
            podman run -e OMP_NUM_THREADS=16 -v /net/co2/c2sm-data/extpar-input-data:/data extpar:$ghprbPullId bash -c "/workspace/test/jenkins/test_docker.sh ${MODE}" || (podman image rm -f extpar:$ghprbPullId && exit 1)
        fi

        podman image rm -f extpar:$ghprbPullId
        exit 0
        ;;
esac
