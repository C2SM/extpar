#!/bin/bash
case "$(hostname)" in
    *co2* | *iacdipl-7*)
        set -e
        podman run -e OMP_NUM_THREADS=16 -v /net/co2/c2sm-data/extpar-input-data:/data extpar:$ghprbPullId bash -c /workspace/test/jenkins/test_docker.sh || (podman image rm -f extpar:$ghprbPullId && exit 1)
        podman image rm -f extpar:$ghprbPullId
        exit 0
        ;;
esac
