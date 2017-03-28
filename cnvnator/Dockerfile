FROM ubuntu:14.04
LABEL maintainer "Colby Chiang <colbychiang@wustl.edu>"

# Build dependencies
RUN apt-get update -qq \
    && apt-get -y install \
        apt-transport-https \
	\
	wget \
	libxpm4 \
	libgomp1 \
	\
       	git \
       	gawk \
        ca-certificates \
    && echo "deb [trusted=yes] https://gitlab.com/hall-lab/ccdg-apt-repo/raw/master ccdg main" | tee -a /etc/apt/sources.list \
    && runDeps=' \
	ccdg-python-2.7.12 \
        ccdg-samtools-1.3.1 \
        ' \
    && apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        $runDeps \
    && rm -rf /var/lib/apt/lists/* \
    && cd /opt \
    && wget -q --no-check-certificate http://colbychiang.com/hall/ccdg-cnvnator-0.3.3_0.3.3-1ubuntu14.04.deb \
    && dpkg --install ccdg-cnvnator-0.3.3_0.3.3-1ubuntu14.04.deb \
    && rm ccdg-cnvnator-0.3.3_0.3.3-1ubuntu14.04.deb

ENV PATH /opt/ccdg/samtools-1.3.1/bin:${PATH}
ENV PATH /opt/ccdg/python-2.7.12/bin:${PATH}

RUN cd /opt \
    && git clone https://github.com/hall-lab/speedseq.git \
    && cd /opt/speedseq \
    && git checkout v0.1.2_cram_support

    # && git submodule sync \
    # && git submodule update --init \
    #     src/bamkit

ENV SHELL /bin/bash

CMD ["/bin/bash"]
