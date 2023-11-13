FROM rocker/r-ver:4.0.3
MAINTAINER Anup Jonchhe <cmap-soft@broadinstitute.org>
LABEL base.mts.pipeline.clue.io.version="0.0.1"
LABEL base.mts.pipeline.clue.io.vendor="PRISM"

RUN apt-get update -qq && \
  apt-get -y upgrade
RUN mkdir -p /src && apt-get -y --no-install-recommends install \
  libssl-dev \
  libxml2-dev \
  libcairo2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadbclient-dev \
  libpq-dev \
  libssh2-1-dev \
  libhdf5-dev \
  libcurl4-openssl-dev \
  jq

RUN apt-get install -y --no-install-recommends libbz2-dev liblzma-dev
COPY ./install_packages.R /src/install_packages.R
RUN mkdir -p /clue/bin && Rscript /src/install_packages.R
COPY ./install_sushi_packages.R /src/install_sushi_packages.R
RUN Rscript /src/install_sushi_packages.R
RUN apt-get -y update; apt-get -y install curl
RUN curl -L -o taigr.zip https://github.com/broadinstitute/taigr/archive/master.zip && \
    unzip taigr.zip && \
    cd taigr-master && \
    R CMD INSTALL .
COPY ./install_sushi_packages2.R /src/install_sushi_packages2.R
RUN Rscript /src/install_sushi_packages2.R
