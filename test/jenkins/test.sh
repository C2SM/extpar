#!/bin/bash
case "$(hostname)" in
    *co2* | *iacdipl-7*)
        set -e

        trap "podman image rm -f extpar:$ghprbPullId || true" EXIT

        MODE=${1:-'NORMAL'}
        if [[ ${MODE} == "DEBUG" ]]; then
            CURRENT_DATETIME=$(date --iso-8601=seconds)
            HASH=$(echo -n "${CURRENT_DATETIME}" | sha256sum | awk '{print $1}')
            mkdir -p /net/co2/c2sm-services/extpar/${HASH}

            if podman run \
               -e OMP_NUM_THREADS=16 \
               -v /net/co2/c2sm-data/extpar-input-data:/data \
               -v /net/co2/c2sm-services/extpar/${HASH}:/artifacts \
               extpar:$ghprbPullId \
               bash -c "/workspace/test/jenkins/test_docker.sh ${MODE}"; then
                REPORT_STATUS="completed successfully"
                STATUS=0
            else
                REPORT_STATUS="failed"
                STATUS=1
            fi

            curl -H "Authorization: token ${GITHUB_TOKEN}" \
                 -H "Accept: application/vnd.github+json" \
                 -H "Content-Type: application/json" \
                 https://api.github.com/repos/C2SM/extpar/issues/${ghprbPullId}/comments \
                 -d "{\"body\": \"The testsuite you launched for this PR has ${REPORT_STATUS}. You can inspect the results at this [link](https://data.iac.ethz.ch/extpar/${HASH}).\"}"
        else
            podman run -e OMP_NUM_THREADS=16 -v /net/co2/c2sm-data/extpar-input-data:/data extpar:$ghprbPullId bash -c "/workspace/test/jenkins/test_docker.sh ${MODE}"
            STATUS="$?"
        fi

        exit "${STATUS}"
        ;;
esac
