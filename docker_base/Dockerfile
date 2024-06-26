FROM rocker/r-ver:4.1.0
MAINTAINER John Davis <cmap-soft@broadinstitute.org>

# Installations

# Set root user for installs
USER root

# Install packages
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

# First set of R dependencies
COPY ./install_packages.R /src/install_packages.R
RUN mkdir -p /clue/bin && Rscript /src/install_packages.R

# Second set of R dependencies
COPY ./install_sushi_packages.R /src/install_sushi_packages.R
RUN Rscript /src/install_sushi_packages.R

# Third set of R dependencies
COPY ./install_sushi_packages2.R /src/install_sushi_packages2.R
RUN Rscript /src/install_sushi_packages2.R

# Install curl (why post update?)
RUN apt-get -y update; apt-get -y install curl

# drop back to the regular jenkins user - good practice
USER jenkins
