FROM debian:stretch-slim
LABEL maintainer "Dave Larson <delarson@wustl.edu>"

COPY --from=halllab/bcftools-1.9-build:v1 /build/deb-build/opt/hall-lab/bcftools-1.9 /opt/hall-lab/bcftools-1.9
COPY --from=halllab/htslib-1.9-build:v1 /build/deb-build/opt/hall-lab/htslib-1.9 /opt/hall-lab/htslib-1.9

ENV PATH=/opt/hall-lab/bcftools-1.9/bin:/opt/hall-lab/htslib-1.9/bin:$PATH

RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        libssl1.1 \
        libcurl3 \
        libbz2-1.0 \ 
        liblzma5 \ 
        libssl1.0.2 \
        zlib1g

CMD ["/bin/bash"]
