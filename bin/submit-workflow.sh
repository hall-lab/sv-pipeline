#!/bin/bash

set -eo pipefail

# current directory of this script
SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
DATA_DIR=${SOURCE_DIR}/..


WDL=${DATA_DIR}/JointGenotyping.wdl
INPUTS=${DATA_DIR}/JointGenotyping.hg38.wgs.inputs.json
OPTIONS=${DATA_DIR}/2020.02.25.options.json
ZIP=${DATA_DIR}/tasks.zip

if [ -e ${ZIP} ]; then
    echo "Nuking out old dependency zip: ${ZIP}"
    rm ${ZIP}
fi

echo "Creating dependency zip: ${ZIP}"
cd ${DATA_DIR}; /usr/bin/zip -r ${ZIP} ./tasks

#echo "WDL=${WDL}"
#echo "INPUTS=${INPUTS}"
#echo "OPTIONS=${OPTIONS}"
#echo "ZIP=${ZIP}"

set -o xtrace

curl -v "localhost:8000/api/workflows/v1" \
	-F workflowSource=@${WDL} \
	-F workflowInputs=@${INPUTS} \
	-F workflowOptions=@${OPTIONS} \
	-F workflowDependencies=@${ZIP}

set +o xtrace

echo "Nuking out dependency zip: ${ZIP}"
rm ${ZIP}

# bash bin/submit-workflow.sh 2>&1 | tee logs/client/submit-job.log.$(date "+%Y.%m.%d.%H.%M")
