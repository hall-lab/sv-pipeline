FROM ubuntu:14.04
LABEL maintainer "Colby Chiang <colbychiang@wustl.edu>"

# Build dependencies
RUN apt-get update -qq \
    && apt-get -y install \
        apt-transport-https \
    && echo "deb [trusted=yes] https://gitlab.com/hall-lab/ccdg-apt-repo/raw/master ccdg main" | tee -a /etc/apt/sources.list \
    && runDeps=' \
	ccdg-python-2.7.12 \
	ccdg-samtools-1.3.1 \
        ' \
    && apt-get update -qq \
    && apt-get -y install \
        --no-install-recommends \
        $runDeps \
    && /opt/ccdg/python-2.7.12/bin/pip install --upgrade pip svtools \
    && rm -rf /var/lib/apt/lists/*

ENV PATH /opt/ccdg/python-2.7.12/bin:${PATH}
ENV PATH /opt/ccdg/samtools-1.3.1/bin:${PATH}

ENV SHELL /bin/bash

CMD ["/bin/bash"]
