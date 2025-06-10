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
    libsqlite3-dev \
    libmariadbd-dev \
    libmariadb-dev-compat \
    libpq-dev \
    libssh2-1-dev \
    libhdf5-dev \
    libcurl4-openssl-dev \
    jq \
    python3 \
    python3-pip && \
    apt-get -y --no-install-recommends install libbz2-dev liblzma-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Python dependencies
RUN pip3 install pandas requests boto3 polars --break-system-packages

# First set of R dependencies
COPY docker/install_packages.R /src/install_packages.R

RUN Rscript /src/install_packages.R

# Install curl (why post update?)
RUN apt-get update -qq && \
    apt-get -y install --no-install-recommends curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
