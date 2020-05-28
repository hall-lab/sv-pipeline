#!/bin/bash

set -eo pipefail

# https://stackoverflow.com/questions/9893667
trap "exit 1" TERM
export TOP_PID=$$

# current directory of this script
SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
DATA_DIR=${SOURCE_DIR}/..

function die {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "[ ${timestamp} ] ERROR: $@" >&2
    kill -s TERM ${TOP_PID}
}

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

function get_input_json {
    local __cohort=$1
    local input_json=${DATA_DIR}/inputs/${__cohort}/${__cohort}.input.json

    if [ ! -e ${input_json} ]; then
        die "Did not find input json file: ${input_json}"
    fi

    echo ${input_json}
}

function get_option_json {
    local __cohort=$1
    local input_option=${DATA_DIR}/inputs/${__cohort}/${__cohort}.options.json

    if [ ! -e ${input_option} ]; then
        die "Did not find option file: ${input_option}"
    fi

    echo ${input_option}
}

function submit_pre_merge_workflow {
    local options=$1
    local inputs=$2
    local scriptsdir=${DATA_DIR}/scripts
    local tmpdir=${DATA_DIR}/tmp
    local wdl=${scriptsdir}/Pre_Merge_SV.wdl
    local zip=tasks.zip
    
    if [ -e ${zip} ]; then
        log "Nuking out old dependency zip: ${zip}"
        rm ${zip}
    fi

    if [ -d ${tmpdir} ]; then
        log "Nuking out tmp workspace: ${tmpdir}"
        rm -rf ${tmpdir}
    fi
    
    log "Creating tmp workspace: ${tmpdir}"
    mkdir -p ${tmpdir}
    cp -v ${scriptsdir}/Pre_Merge_SV.wdl ${tmpdir}/
    cp -v ${scriptsdir}/Pre_Merge_SV_per_sample.wdl ${tmpdir}/
    cp -v ${scriptsdir}/Pre_Merge_QC_per_sample.wdl ${tmpdir}/
    cp -v ${scriptsdir}/SV_Tasks.wdl ${tmpdir}/
    log "Creating dependency zip: ${zip}"
    cd ${tmpdir} && /usr/bin/zip -r tasks.zip Pre_Merge_QC_per_sample.wdl Pre_Merge_SV_per_sample.wdl SV_Tasks.wdl
    cd ${DATA_DIR}; mv ${tmpdir}/tasks.zip ${DATA_DIR}
    
    #log "WDL=${wdl}"
    #log "INPUTS=${inputs}"
    #log "OPTIONS=${options}"
    #log "ZIP=${zip}"
    
    set -o xtrace
    
    curl -v "localhost:8000/api/workflows/v1" \
    	-F workflowSource=@${wdl} \
    	-F workflowInputs=@${inputs} \
    	-F workflowOptions=@${options} \
    	-F workflowDependencies=@${zip}
    
    set +o xtrace
    
    log "Nuking out dependency zip: ${zip}"
    rm ${zip}
    log "Nuking out tmp workspace: ${tmpdir}"
    rm -rf ${tmpdir}
}

function main {
    local cohort=$1
    local option_json=$(get_option_json ${cohort})
    local input_json=$(get_input_json ${cohort})
    submit_pre_merge_workflow ${option_json} ${input_json}
}

COHORT=$1;
main ${COHORT};

## mkdir -p logs/client # if the logs/client direct doesn't exist in the repo workspace yet
## COHORT=afib
## bash bin/submit-workflow.sh ${COHORT} 2>&1 | tee logs/client/submit-job.log.${COHORT}.$(date "+%Y.%m.%d.%H.%M")
