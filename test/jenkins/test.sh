#!/bin/bash
case "$(hostname)" in
    *co2* | *iacdipl-7*)
        set -e
        podman run -e OMP_NUM_THREADS=16 -v /net/co2/c2sm-data/extpar-input-data:/data -v /net/co2/c2sm-services/extpar/test:/artifacts extpar:$ghprbPullId bash -c "/workspace/test/jenkins/test_docker.sh ${1:-'2'}" || (podman image rm -f extpar:$ghprbPullId && exit 1)
        podman image rm -f extpar:$ghprbPullId
        exit 0
        ;;
esac
