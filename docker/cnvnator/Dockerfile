FROM debian:stretch-slim
LABEL maintainer "Dave Larson <delarson@wustl.edu>"
COPY --from=halllab/python2.7-build:v1 /opt/hall-lab/python-2.7.15 /opt/hall-lab/python-2.7.15
COPY --from=halllab/htslib-1.9-build:v1 /build/deb-build/opt/hall-lab/htslib-1.9 /opt/hall-lab/htslib-1.9
COPY --from=halllab/samtools-1.9-build:v1 /build/deb-build/opt/hall-lab/samtools-1.9 /opt/hall-lab/samtools-1.9
COPY --from=halllab/cnvnator-0.3.3-build:v1 /opt/hall-lab/cnvnator-0.3.3/deb-build/opt/hall-lab/cnvnator-0.3.3 /opt/hall-lab/cnvnator-0.3.3
ENV PATH=/opt/hall-lab/python-2.7.15/bin:/opt/hall-lab/htslib-1.9/bin:/opt/hall-lab/samtools-1.9/bin:/opt/hall-lab/cnvnator-0.3.3/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/hall-lab/htslib-1.9/lib:$LD_LIBRARY_PATH

# Build dependencies
RUN apt-get update -qq \
    && apt-get -y install \
        libssl1.1 \
        libcurl3 \
        libncurses5 \
        libbz2-1.0 \ 
        liblzma5 \ 
        libssl1.0.2 \
        zlib1g
        libgomp1 \
        libstdc++6 \
        libstdc++-6-dev \
        libgcc1 \
        g++ \
        libxpm4 \
        git-core \
    && git clone https://github.com/hall-lab/speedseq.git \
    && cd speedseq \
    && git checkout 4e60002 \
    && cp bin/cnvnator_wrapper.py /opt/hall-lab/cnvnator-0.3.3/bin/ \
    && cd .. \
    && rm -rf speedseq \
    && apt-get purge -y git-core \
    && apt autoremove -y

CMD ["/bin/bash"]
