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
        ca-certificates \
        tabix \
        libbz2-dev \
        liblzma-dev \
    && pip install pandas \
    && pip install scipy \
    && pip install pysam \
    && pip install svtools
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
LABEL description "Manta v1.4.0 with samtools v1.9 alongside"
ARG MANTA_VERSION=1.4.0

COPY --from=manta-build /opt/hall-lab/manta-${MANTA_VERSION}/bin /opt/hall-lab/manta-${MANTA_VERSION}/bin
COPY --from=manta-build /opt/hall-lab/manta-${MANTA_VERSION}/lib /opt/hall-lab/manta-${MANTA_VERSION}/lib
COPY --from=manta-build /opt/hall-lab/manta-${MANTA_VERSION}/libexec /opt/hall-lab/manta-${MANTA_VERSION}/libexec
COPY --from=manta-build /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
COPY --from=halllab/htslib-1.9-build:v1 /build/deb-build/opt/hall-lab/htslib-1.9 /opt/hall-lab/htslib-1.9
COPY --from=halllab/samtools-1.9-build:v1 /build/deb-build/opt/hall-lab/samtools-1.9 /opt/hall-lab/samtools-1.9
ADD doctor_manta.1.py /opt/hall-lab/

# Run dependencies
RUN apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        libssl1.1 \
        libcurl3 \
        libncurses5 \
        libbz2-1.0 \ 
        liblzma5 \ 
        libssl1.0.2 \
        zlib1g

ENV PATH=/opt/hall-lab/manta-${MANTA_VERSION}/bin:/opt/hall-lab/python-2.7.15/bin/:/opt/hall-lab/samtools-1.9/bin:/opt/hall-lab/htslib-1.9/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/hall-lab/htslib-1.9/lib:$LD_LIBRARY_PATH

CMD ["/bin/bash"]
