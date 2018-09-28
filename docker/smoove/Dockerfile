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


FROM builder-base as lumpy-ea9691f-build
LABEL maintainer "Dave Larson <delarson@wustl.edu>"
RUN LUMPY_COMMIT=ea9691f52adcb38cfad8edc1f29a058c581297c6 \
    && git clone --single-branch --recursive --depth 5 https://github.com/arq5x/lumpy-sv \
    && cd lumpy-sv \
    && git checkout $LUMPY_COMMIT \
    && make -j 3 \
    && mkdir -p /opt/hall-lab/lumpy-ea9691f/bin \
    && cp ./bin/* /opt/hall-lab/lumpy-ea9691f/bin
    
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
RUN find /opt/hall-lab/python-2.7.15/lib/python2.7/site-packages/ -name '*.so' -print -exec sh -c 'file "{}" | grep -q "not stripped" && strip -s "{}"' \;

FROM builder-base as extract-sv-reads-build
RUN LIBDEFLATE_VERSION=1.0 \
    && EXTRACT_SV_VERSION=1.3.0 \
    && git clone https://github.com/ebiggers/libdeflate.git \
    && cd libdeflate \
    && git checkout v${LIBDEFLATE_VERSION} \
    && make -j 2 CFLAGS='-fPIC -O3' libdeflate.a \
    && install -Dm644 -t /usr/local/lib libdeflate.a \
    && install -Dm644 -t /usr/local/include libdeflate.h \
    && ldconfig \
    && cd .. \
    && git clone -b v${EXTRACT_SV_VERSION} https://github.com/hall-lab/extract_sv_reads \
    && mkdir build \
    && cd build \
    && cmake -DHTSLIB_USE_LIBCURL=1 -DHTSLIB_USE_LZMA=1 -DHTSLIB_USE_BZ2=1 -DHTSLIB_USE_LIBDEFLATE=1 $WORKDIR/extract_sv_reads/ \
    && make -j 2 \
    && mkdir -p /opt/hall-lab/extract-sv-reads-${EXTRACT_SV_VERSION}/bin \
    && cp bin/extract-sv-reads /opt/hall-lab/extract-sv-reads-${EXTRACT_SV_VERSION}/bin \
    && git clone https://gist.github.com/ernfrid/9cf947d9f9d2582be291399862c96b86 \
    && cp 9cf947d9f9d2582be291399862c96b86/lumpy_filter /opt/hall-lab/extract-sv-reads-${EXTRACT_SV_VERSION}/bin/ \
    && chmod a+x /opt/hall-lab/extract-sv-reads-${EXTRACT_SV_VERSION}/bin/lumpy_filter

# Smoove build...
FROM builder-base as smoove-0.1.10-build
WORKDIR /opt/hall-lab/smoove-0.1.10/bin
RUN SMOOVE_VERSION=0.1.10 \
    && MOSDEPTH_VERSION=0.2.3 \
    && GSORT_VERSION=0.0.6 \
    && curl -L -o mosdepth https://github.com/brentp/mosdepth/releases/download/v$MOSDEPTH_VERSION/mosdepth \
    && chmod a+x mosdepth \
    && curl -L -o gsort https://github.com/brentp/gsort/releases/download/v$GSORT_VERSION/gsort_linux_amd64 \
    && chmod a+x gsort \
    && curl -L -o smoove https://github.com/brentp/smoove/releases/download/v$SMOOVE_VERSION/smoove \
    && chmod +x smoove


FROM debian:stretch-slim
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

COPY --from=lumpy-ea9691f-build /opt/hall-lab/lumpy-ea9691f/bin /opt/hall-lab/lumpy-ea9691f/bin
COPY --from=svtyper-0.7.0-build /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
COPY --from=halllab/htslib-1.9-build:v1 /build/deb-build/opt/hall-lab/htslib-1.9 /opt/hall-lab/htslib-1.9
COPY --from=halllab/samtools-1.9-build:v1 /build/deb-build/opt/hall-lab/samtools-1.9 /opt/hall-lab/samtools-1.9
COPY --from=halllab/bcftools-1.9-build:v1 /build/deb-build/opt/hall-lab/bcftools-1.9 /opt/hall-lab/bcftools-1.9
COPY --from=smoove-0.1.10-build /opt/hall-lab/smoove-0.1.10/bin /opt/hall-lab/smoove-0.1.10/bin
COPY --from=extract-sv-reads-build /opt/hall-lab/extract-sv-reads-1.3.0/bin /opt/hall-lab/extract-sv-reads-1.3.0/bin

ENV PATH=/opt/hall-lab/smoove-0.1.10/bin:/opt/hall-lab/python-2.7.15/bin:/opt/hall-lab/extract-sv-reads-1.3.0/bin:/opt/hall-lab/lumpy-ea9691f/bin:/opt/hall-lab/htslib-1.9/bin:/opt/hall-lab/samtools-1.9/bin:/opt/hall-lab/bcftools-1.9/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/hall-lab/htslib-1.9/lib:$LD_LIBRARY_PATH

RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        libssl1.1 \
        libcurl3 \
        libncurses5 \
        libbz2-1.0 \ 
        liblzma5 \ 
        libssl1.0.2 \
        zlib1g

CMD ["/bin/bash"]
