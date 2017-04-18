JES_CONF=jes.conf
CROMWELL=/home/cchiang/src/cromwell/target/scala-2.11/cromwell-26-22fe860-SNAP.jar

java \
    -Dconfig.file=$JES_CONF \
    -jar $CROMWELL \
    run \
    test.wdl \
    test.inputs.json \
    ../../scripts/generic.options.json \
    test.metadata.json
