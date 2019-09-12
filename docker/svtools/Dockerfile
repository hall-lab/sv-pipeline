FROM debian:stretch-slim as builder-base
LABEL maintainer "Dave Larson <delarson@wustl.edu>"
RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        cmake \
        autoconf \
        automake \
        libtool \
        gawk \
        git-core \
        bzip2 \
        libbz2-dev \
        liblzma-dev \
        libssl1.0-dev \
        libcurl4-openssl-dev \
        ca-certificates \
        libblas-dev \
        libatlas-base-dev \
        liblapack-dev \
        curl \
        zlib1g-dev

FROM builder-base as svtools-0.5.1-build
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

COPY --from=halllab/python2.7-build:v1 /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
ENV PATH=/opt/hall-lab/python-2.7.15/bin:${PATH}
RUN SVTOOLS_VERSION=0.5.1 \
    && pip install svtools==${SVTOOLS_VERSION}
RUN find /opt/hall-lab/python-2.7.15/ -depth \( -name '*.pyo' -o -name '*.pyc' -o -name 'test' -o -name 'tests' \) -exec rm -rf '{}' + ;
#RUN find /opt/hall-lab/python-2.7.15/lib/python2.7/site-packages/ -name '*.so' -print -exec sh -c 'file "{}" | grep -q "not stripped" && strip -s "{}"' \;

FROM debian:stretch-slim
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

COPY --from=svtools-0.5.1-build /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
COPY --from=halllab/samtools-1.9-build:v1 /build/deb-build/opt/hall-lab/samtools-1.9 /opt/hall-lab/samtools-1.9
COPY --from=halllab/cnvnator-0.3.3-build:v1 /opt/hall-lab/cnvnator-0.3.3/deb-build/opt/hall-lab/cnvnator-0.3.3 /opt/hall-lab/cnvnator-0.3.3
COPY --from=halllab/htslib-1.9-build:v1 /build/deb-build/opt/hall-lab/htslib-1.9 /opt/hall-lab/htslib-1.9

RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        libssl1.1 \
        libcurl3 \
        libncurses5 \
        libbz2-1.0 \ 
        liblzma5 \ 
        libssl1.0.2 \
        zlib1g \
        libblas3 \
        libatlas3-base \
        liblapack3 \
        libgomp1 \
        libstdc++6 \
        libstdc++-6-dev \
        libgcc1 \
        g++ \
        libxpm4 \
        gzip \
        less

ENV PATH=/opt/hall-lab/python-2.7.15/bin:/opt/hall-lab/samtools-1.9/bin:/opt/hall-lab/htslib-1.9/bin:/opt/hall-lab/cnvnator-0.3.3/bin:${PATH}

CMD ["/bin/bash"]

