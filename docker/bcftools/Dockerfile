FROM ubuntu:14.04
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

# Build dependencies
RUN apt-get update -qq \
    && apt-get -y install apt-transport-https \
    && echo "deb [trusted=yes] https://gitlab.com/hall-lab/ccdg-apt-repo/raw/master ccdg main" | tee -a /etc/apt/sources.list \
    && runDeps=' \
        libcurl3 \
        ca-certificates \
        zlib1g \
        libncurses5 \
        ccdg-bcftools-1.3.1 \
        ' \
    && apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        $runDeps \
    && rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/ccdg/bcftools-1.3.1/bin:${PATH}

CMD ["/bin/bash"]
