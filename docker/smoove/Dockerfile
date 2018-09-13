FROM ubuntu:16.04
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

# Build dependencies
RUN SMOOVE_VERSION=0.1.7 \
    && SAMTOOLS_VERSION=1.8 \
    && BCFTOOLS_VERSION=1.8 \
    && LIBDEFLATE_VERSION=1.0 \
    && SVTYPER_VERSION=0.7.0 \
    && MOSDEPTH_VERSION=0.2.3 \
    && GSORT_VERSION=0.0.6 \
    && LUMPY_COMMIT=ea9691f52adcb38cfad8edc1f29a058c581297c6 \
    && buildDeps=' \
        build-essential \
        make \
        autoconf \
        automake \
        libtool \
        gawk \
        git-core \
        bzip2 \
        libbz2-dev \
        liblzma-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libffi-dev \
        libblas-dev \
        liblapack-dev \
        libatlas-base-dev \
        libuuid1 \
        curl \
        ca-certificates \
        zlib1g \
        zlib1g-dev \
        libncurses-dev \
        python \
        python-dev \
        python-pip \
        ' \
    && runDeps=' \
        python \
        zlib1g \
        ' \
    && apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        $buildDeps \
    && git clone https://github.com/ebiggers/libdeflate.git \
    && cd libdeflate \
    && git checkout v$LIBDEFLATE_VERSION \
    && make -j 2 CFLAGS='-fPIC -O3' libdeflate.a \
    && cp libdeflate.a /usr/local/lib \
    && cp libdeflate.h /usr/local/include \
    && cd .. && rm -rf libdeflate \
    && curl -OL https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 \
    && tar -xjf samtools-$SAMTOOLS_VERSION.tar.bz2 \
    && cd samtools-$SAMTOOLS_VERSION \
    && ./configure --enable-libcurl --with-libdeflate \
    && make -j2 all \
    && make -j2 install \
    && cd htslib-$SAMTOOLS_VERSION && make install \
    && cd .. \
    && cd .. && rm -rf samtools-$BCFTOOLS_VERSION \
    && curl -OL https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/bcftools-$BCFTOOLS_VERSION.tar.bz2 \
    && tar -xjf bcftools-$BCFTOOLS_VERSION.tar.bz2 \
    && cd bcftools-$BCFTOOLS_VERSION \
    && ./configure --enable-libcurl --with-libdeflate \
    && make -j2 all \
    && make -j2 install \
    && cd .. && rm -rf bcftools-$BCFTOOLS_VERSION \
    && pip install setuptools \
    && git clone https://github.com/hall-lab/svtyper \
    && cd svtyper \
    && git checkout v$SVTYPER_VERSION \
    && python setup.py install \
    && cd .. && rm -rf svtyper \
    && curl -L -o /usr/local/bin/mosdepth https://github.com/brentp/mosdepth/releases/download/v$MOSDEPTH_VERSION/mosdepth \
    && chmod a+x /usr/local/bin/mosdepth \
    && curl -L -o /usr/bin/gsort https://github.com/brentp/gsort/releases/download/v$GSORT_VERSION/gsort_linux_amd64 \
    && chmod a+x /usr/bin/gsort \
    && git clone --single-branch --recursive --depth 1 https://github.com/arq5x/lumpy-sv \
    && cd lumpy-sv \
    && git checkout $LUMPY_COMMIT \
    && make -j 3 \
    && cp ./bin/* /usr/local/bin/ \
    && cd .. && rm -rf lumpy-sv \
    && curl -L -o /usr/local/bin/smoove https://github.com/brentp/smoove/releases/download/v$SMOOVE_VERSION/smoove \
    && chmod +x /usr/local/bin/smoove \
    && ldconfig \
    && rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]
