JES_CONF=jes.conf
CROMWELL=/home/cchiang/src/cromwell/target/scala-2.11/cromwell-26-22fe860-SNAP.jar
OPTIONS=../../../scripts/generic.options.json

java \
    -Dconfig.file=$JES_CONF \
    -jar $CROMWELL \
    run \
    test.wdl \
    test.inputs.json \
    $OPTIONS \
    test.metadata.json
