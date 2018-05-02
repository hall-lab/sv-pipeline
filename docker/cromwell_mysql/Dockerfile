FROM broadinstitute/cromwell:29
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

# Build dependencies
RUN apt-get update -qq \
    && runDeps=' \
        libnss-sss \
        mysql-server \
        ' \
    && apt-get update -qq \
    && DEBIAN_FRONTEND=noninteractive apt-get -y install \
        --no-install-recommends \
        $runDeps \
    && mkdir /var/run/mysqld \
    && chmod -R 777 /var/lib/mysql /var/run/mysqld /var/log/mysql && rm -fr /var/lib/mysql/mysql /var/lib/mysql/performance_schema \
    && mkdir -p /opt/ccdg/cromwell/resources \
    && ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime && echo "America/Chicago" > /etc/timezone && dpkg-reconfigure --frontend noninteractive tzdata \
    && rm -rf /var/lib/apt/lists/*

ADD application.conf.template /opt/ccdg/cromwell/resources
ADD mysql.cnf.template /opt/ccdg/cromwell/resources
ADD run_pipeline.sh /opt/ccdg/cromwell/resources

# Reset entrypoint so container doesn't try to run cromwell directly
ENTRYPOINT ["/opt/ccdg/cromwell/resources/run_pipeline.sh"]
