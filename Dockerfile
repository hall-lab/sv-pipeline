FROM ubuntu:14.04
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

# Build dependencies
RUN export EXTRACT_SV_READS_VERSION=1.1.0 \
    && apt-get update -qq \
    && apt-get -y install apt-transport-https \
    && echo "deb [trusted=yes] https://gitlab.com/hall-lab/ccdg-apt-repo/raw/master ccdg main" | tee -a /etc/apt/sources.list \
    && runDeps=' \
        libcurl3 \
        ca-certificates \
        zlib1g \
        libncurses5 \
        ccdg-samtools-1.3.1 \
        extract-sv-reads1.1 \
        ' \
    && apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        $runDeps \
    && rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/ccdg/samtools-1.3.1/bin:${PATH}

CMD ["/bin/bash"]
