# Cromwell Mechanics

## Cromwell Server Setup

The cromwell server process is done via the [cromwell-deployment][0] repository.  The particular custom configurations for the sv-pipeline are placed in the [config/cromwell-server][1].  Please see the [decryption instructions][2] on how to obtain the custom configuration file.

### Step 1:  Setup `cromwell-deployment` &amp; Invoke the Google Deployment

    # assuming that you are already in the root of this github repository
    cd config/cromwell-server
    
    # decrypt the configuration file
    gpg --output cromwell.sv.yaml --no-symkey-cache --decrypt cromwell.sv.yaml.gpg
    
    # checkout the cromwell-deployment repository
    git clone https://github.com/hall-lab/cromwell-deployment cromwell-deployment
    
    # setup the cromwell-deployment system
    cd cromwell-deployment/resources/google-deployments
    cp ../../../cromwell.sv.yaml .
    
    # invoke the google deployment manager
    gcloud deployment-manager deployments create cromwell1 --config cromwell.sv.yaml

See [cromwell-deployment][0] for more details and instructions on customizing the cromwell setup process.

### Step 2: Setup the Instance

Once the compute instance is available, in this case `cromwell-sv-20200508`, go ahead an ssh into the box as:

    gcloud compute ssh --zone us-central1-c cromwell-sv-20200508

Run the following commands on the box:

    # the location where we'll install all the custom software needed for the sv-pipeline
    cd /opt
    umask 002                 # allow group write; everyone must do this
    chgrp google-sudoers .    # set directory group to GROUPNAME (google-sudoers)
    chmod g+s .               # files created in directory will be in group GROUPNAME (google-sudoers)

    # install a custom version of python (maintained by hall lab members)
    sudo apt-get -y install apt-transport-https
    echo "deb [trusted=yes] https://gitlab.com/indraniel/hall-lab-debian-repo-1/raw/master stretch main" | sudo tee -a /etc/apt/sources.list
    sudo apt-get update -qq
    sudo apt-get -y install hall-lab-python-3.7.0

    # install some other dependencies needed to run other useful programs
    sudo apt-get -y install zip tmux git jq mailutils bsdmainutils curl libcurl4-openssl-dev gnupg

    # obtain the relevant sv-pipeline code
    git clone https://github.com/hall-lab/sv-pipeline
    cd sv-pipeline
    git checkout -b post-freeze-2-2020-05-08 origin/post-freeze-2-2020-05-08

    # setup a virtualenv for running certain python scripts
    /opt/hall-lab/python-3.7.0/bin/virtualenv venv37
    source venv37/bin/activate
    pip install -r requirements.txt

    # install cromshell
    cd /opt
    git clone https://github.com/broadinstitute/cromshell
    cd cromshell
    ./cromshell -h

## Workflow Submissions

### Step 1:  Generate the json options and inputs

    gcloud compute ssh --zone us-central1-c cromwell-sv-20200508
    cd /opt
    umask 002                 # allow group write; everyone must do this
    sudo chgrp google-sudoers .    # set directory group to GROUPNAME (google-sudoers)
    sudo chmod g+s .               # files created in directory will be in group GROUPNAME (google-sudoers)
    cd /opt/sv-pipeline
    source venv37/bin/activate
    python bin/make-json.py --json-type="input" --sample-map=data/foo.tsv --cohort=afib
    python bin/make-json.py --json-type="option" --cohort=afib
    # you should see the json files inside the input/<cohort> directory of the repository

### Step 2:  Submit the workflow

    gcloud compute ssh --zone us-central1-c cromwell-sv-20200508
    cd /opt
    umask 002                 # allow group write; everyone must do this
    sudo chgrp google-sudoers .    # set directory group to GROUPNAME (google-sudoers)
    sudo chmod g+s .               # files created in directory will be in group GROUPNAME (google-sudoers)
    cd /opt/sv-pipeline
    mkdir -p logs/client # if the logs/client direct doesn't exist in the repo workspace yet
    COHORT=afib
    bash bin/submit-workflow.sh ${COHORT} 2>&1 | tee logs/client/submit-job.log.${USER}.$(date "+%Y.%m.%d.%H.%M")

Afterwards, note the output workflow id.  You'll need the workflow id to monitor status and progress.

## Monitoring

    gcloud compute ssh --zone us-central1-c cromwell-sv-20200508 --command="/opt/sv-pipeline/bin/status-workflow.sh <workflow-id>"
    ## or
    gcloud compute ssh --zone us-central1-c cromwell-sv-20200508 --command="/opt/cromshell/cromshell status <workflow-id>"

To get the detailed status information try:

    glcoud compute ssh --zone us-central1-c cromwell-sv-20200508 --command="/opt/cromshell/cromshell execution-status-count -p -x <workflow-id>"

[0]:  https://github.com/hall-lab/cromwell-deployment
[1]:  https://github.com/hall-lab/sv-pipeline/tree/post-freeze-2-2020-05-08/config/cromwell-server
[2]:  https://github.com/hall-lab/sv-pipeline/blob/post-freeze-2-2020-05-08/config/cromwell-server/encryption-note.md
