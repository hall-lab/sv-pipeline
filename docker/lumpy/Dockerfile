FROM ubuntu:14.04
LABEL maintainer "Colby Chiang <colbychiang@wustl.edu>"

# Build dependencies
RUN apt-get update -qq \
    && apt-get -y install \
        apt-transport-https \
        g++ \
	gawk \
        libcurl4-gnutls-dev \
        autoconf \
	libssl-dev \
        git \
    && echo "deb [trusted=yes] https://gitlab.com/hall-lab/ccdg-apt-repo/raw/master ccdg main" | tee -a /etc/apt/sources.list \
    && runDeps=' \
	ccdg-python-2.7.12 \
        ccdg-samtools-1.3.1 \
        ' \
    && apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        $runDeps \
    && /opt/ccdg/python-2.7.12/bin/pip install --upgrade pip numpy scipy pysam \
    && rm -rf /var/lib/apt/lists/*

ENV PATH /opt/ccdg/samtools-1.3.1/bin:${PATH}
ENV PATH /opt/ccdg/python-2.7.12/bin:${PATH}

# Install LUMPY
RUN cd /opt \
    && git clone https://github.com/hall-lab/lumpy-sv.git \
    && cd /opt/lumpy-sv \
    && git checkout 0.2.13_cram_support \
    && git submodule sync \
    && git submodule update --init \
    && cd /opt/lumpy-sv \
    && make

ENV PATH /opt/lumpy-sv/bin:${PATH}
ENV SHELL /bin/bash

CMD ["/bin/bash"]
