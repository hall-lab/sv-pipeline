FROM ubuntu:14.04
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

# Build dependencies
RUN export EXTRACT_SV_READS_VERSION=1.1.0 \
    && apt-get update -qq \
    && apt-get -y install \
       apt-transport-https \
       \
       wget \
       build-essential \
       zlib1g-dev \
       libncurses5-dev \
       libbz2-dev \
       liblzma-dev \
       git \
       bsdmainutils \
       gawk \
       \
    && echo "deb [trusted=yes] https://gitlab.com/hall-lab/ccdg-apt-repo/raw/master ccdg main" | tee -a /etc/apt/sources.list \
    && runDeps=' \
        libcurl3 \
        ca-certificates \
        zlib1g \
        libncurses5 \
	\
	ccdg-python-2.7.12 \
	\
        ccdg-samtools-1.3.1 \
        extract-sv-reads1.1 \
        ' \
    && apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        $runDeps \
    && rm -rf /var/lib/apt/lists/* \
    \
    \
    && echo "deb [trusted=yes] https://gitlab.com/hall-lab/ccdg-apt-repo/raw/master ccdg main" | tee -a /etc/apt/sources.list \
    && runDeps=' \
       libc6 \
       libgomp1 \
       libstdc++6 \
       libgcc1 \
       libxpm4 \
       ccdg-python-2.7.12 \
       ' \
    && apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        $runDeps \
    && rm -rf /var/lib/apt/lists/* \
    \
    && cd /opt \
    && wget -q --no-check-certificate http://colbychiang.com/hall/ccdg-cnvnator-0.3.3_0.3.3-1ubuntu14.04.deb \
    && dpkg --install ccdg-cnvnator-0.3.3_0.3.3-1ubuntu14.04.deb \
    && rm ccdg-cnvnator-0.3.3_0.3.3-1ubuntu14.04.deb \
    \
    && cd /opt \
    &&  wget -q --no-check-certificate https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 \
    && tar -xf samtools-1.4.tar.bz2 && rm samtools-1.4.tar.bz2 \
    && cd samtools-1.4 \
    && make \
    && make install \
    && make clean

ENV PATH /opt/ccdg/samtools-1.3.1/bin:${PATH}
ENV PATH /opt/ccdg/python-2.7.12/bin:${PATH}

############################################

# speedseq

RUN cd /opt \
    && git clone https://github.com/ernfrid/speedseq.git \
    && cd /opt/speedseq \
    && git checkout v0.1.2_cram_support \
    && git submodule sync \
    && git submodule update --init \
        src/bamkit \
	src/lumpy-sv \
	src/parallel \
	src/samblaster \
	src/svtyper \
	src/tabix \
	src/vawk \
    \
    && cd /opt/speedseq \
    && make sv \
    \
    && cd /opt/speedseq/bin \
    && scp /opt/ccdg/cnvnator-0.3.3/bin/cnvnator . \
    && ln -s /opt/ccdg/cnvnator-0.3.3/bin/cnvnator cnvnator-multi


ENV PATH /opt/speedseq/bin:${PATH}
ENV SHELL /bin/bash

RUN pip install --upgrade pip numpy scipy pysam svtools

ENV REF_PATH /data/.cache/%2s/%2s/%s
ENV REF_CACHE /data/.cache/%2s/%2s/%s

#####################################
# svtools



CMD ["/bin/bash"]
