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
        curl \
        zlib1g-dev

FROM builder-base as svtyper-0.7.0-build
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

COPY --from=halllab/python2.7-build:v1 /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
ENV PATH=/opt/hall-lab/python-2.7.15/bin:${PATH}
RUN SVTYPER_VERSION=0.7.0 \
    && git clone https://github.com/hall-lab/svtyper \
    && cd svtyper \
    && git checkout v$SVTYPER_VERSION \
    && sed -i '/numpy/d' setup.py \
    && sed -i '/scipy/d' setup.py \
    && pip install .
RUN find /opt/hall-lab/python-2.7.15/ -depth \( -name '*.pyo' -o -name '*.pyc' -o -name 'test' -o -name 'tests' \) -exec rm -rf '{}' + ;
#RUN find /opt/hall-lab/python-2.7.15/lib/python2.7/site-packages/ -name '*.so' -print -exec sh -c 'file "{}" | grep -q "not stripped" && strip -s "{}"' \;

FROM debian:stretch-slim
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

COPY --from=svtyper-0.7.0-build /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
COPY --from=halllab/samtools-1.9-build:v1 /build/deb-build/opt/hall-lab/samtools-1.9 /opt/hall-lab/samtools-1.9
COPY --from=halllab/htslib-1.9-build:v1 /build/deb-build/opt/hall-lab/htslib-1.9 /opt/hall-lab/htslib-1.9

RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        libssl1.1 \
        libcurl3 \
        libncurses5 \
        libbz2-1.0 \ 
        liblzma5 \ 
        libssl1.0.2 \
        zlib1g

ENV PATH=/opt/hall-lab/python-2.7.15/bin:/opt/hall-lab/samtools-1.9/bin:/opt/hall-lab/htslib-1.9/bin:${PATH}

CMD ["/bin/bash"]
