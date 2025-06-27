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

# Setup venvironment for Python
RUN apt-get update \
 && apt-get install -y python3-venv python3-pip \
 && python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --upgrade pip setuptools wheel

# install package & deps inside the venv
COPY scripts/utils/prism_tools ./prism_tools
RUN /opt/venv/bin/pip install -e prism_tools \
 && /opt/venv/bin/pip install pandas requests boto3 polars

# Copy the R scripts to the container (for use with AWS)
COPY scripts/drc /scripts/drc
COPY scripts/biomarker /scripts/biomarker
COPY scripts/utils /scripts/utils

# make the venv’s python/pip the default
ENV PATH="/opt/venv/bin:${PATH}"

# First set of R dependencies
COPY docker/install_packages.R /src/install_packages.R

RUN Rscript /src/install_packages.R

# Install curl (why post update?)
RUN apt-get update -qq && \
    apt-get -y install --no-install-recommends curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
