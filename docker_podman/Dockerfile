FROM rocker/r-ver:4.4

MAINTAINER John Davis <cmap-soft@broadinstitute.org>

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Set root user for installs
USER root

# Install required packages
RUN apt-get update -qq && \
    mkdir -p /src && \
    apt-get -y --no-install-recommends install \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libsqlite-dev \
    libmariadbd-dev \
    libmariadb-dev-compat \
    libpq-dev \
    libssh2-1-dev \
    libhdf5-dev \
    libcurl4-openssl-dev \
    jq && \
    apt-get -y --no-install-recommends install libbz2-dev liblzma-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# First set of R dependencies
COPY ./install_packages.R /src/install_packages.R
RUN Rscript /src/install_packages.R

# Install curl (why post update?)
RUN apt-get update -qq && \
    apt-get -y install --no-install-recommends curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

