# samtools

# build:  docker build -t skyport/samtools:1.2 .
# run:    docker run --rm -t -i  skyport/samtools:1.2 samtools --help


FROM ubuntu:14.04

RUN apt-get update && apt-get -y install \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev

ADD https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 /

RUN tar xvjf samtools-1.4.tar.bz2 && \
    cd /samtools-1.4 && \
    make && \
    mv /samtools-1.4/samtools /usr/local/bin/

#ENTRYPOINT ["/usr/local/bin/samtools"]
#CMD ["--help"]
