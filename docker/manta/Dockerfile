FROM debian:stretch-slim AS manta-build
LABEL maintainer "Dave Larson <delarson@wustl.edu>"
ARG MANTA_VERSION=1.4.0
COPY --from=halllab/python2.7-build:v1 /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
ENV PATH=/opt/hall-lab/python-2.7.15/bin:${PATH}
RUN apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        build-essential \
        bzip2 \
        zlib1g-dev \
        curl \
        ca-certificates
RUN curl -O -L https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2 \
    && tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2 \
    && mkdir build \
    && cd build \
    && ../manta-${MANTA_VERSION}.release_src/configure --prefix=/opt/hall-lab/manta-${MANTA_VERSION} \
    && make -j 4 install
RUN find /opt/hall-lab/python-2.7.15/ -depth \( -name '*.pyo' -o -name '*.pyc' -o -name 'test' -o -name 'tests' \) -exec rm -rf '{}' + ;
RUN find /opt/hall-lab/python-2.7.15/lib/python2.7/site-packages/ -name '*.so' -print -exec sh -c 'file "{}" | grep -q "not stripped" && strip -s "{}"' \;

FROM debian:stretch-slim
LABEL maintainer "Dave Larson <delarson@wustl.edu>"
ARG MANTA_VERSION=1.4.0

COPY --from=manta-build /opt/hall-lab/manta-${MANTA_VERSION}/[^s]* /opt/hall-lab/manta-${MANTA_VERSION}/
COPY --from=manta-build /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15

# Run dependencies
RUN apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        zlib1g

ENV PATH=/opt/hall-lab/manta-${MANTA_VERSION}/bin:/opt/hall-lab/python-2.7.15/bin/:$PATH

CMD ["/bin/bash"]
