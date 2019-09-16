FROM debian:stretch-slim as builder-base
LABEL maintainer "Allison Regier <aregier@wustl.edu>"
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

FROM builder-base as vawk-build
LABEL maintainer "Allison Regier <aregier@wustl.edu>"

COPY --from=halllab/python2.7-build:v1 /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
ENV PATH=/opt/hall-lab/python-2.7.15/bin:${PATH}
RUN git clone https://github.com/cc2qe/vawk
RUN git clone https://github.com/hall-lab/io
RUN curl -kL https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools > bedtools
RUN chmod a+x bedtools

FROM debian:stretch-slim
LABEL maintainer "Allison Regier <aregier@wustl.edu>"

COPY --from=vawk-build /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
COPY --from=halllab/htslib-1.9-build:v1 /build/deb-build/opt/hall-lab/htslib-1.9 /opt/hall-lab/htslib-1.9
COPY --from=vawk-build vawk /opt/hall-lab/vawk
COPY --from=vawk-build io /opt/hall-lab/io
COPY --from=vawk-build bedtools /opt/hall-lab/bin/bedtools

RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        libssl1.1 \
        libcurl3 \
        libncurses5 \
        libbz2-1.0 \ 
        liblzma5 \ 
        libssl1.0.2 \
        zlib1g \
        less \
	gawk

ENV PATH=/opt/hall-lab/python-2.7.15/bin:${PATH}

CMD ["/bin/bash"]
