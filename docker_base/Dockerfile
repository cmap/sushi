FROM rocker/r-ver:4.1.0
MAINTAINER John Davis <cmap-soft@broadinstitute.org>

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Set root user for installs
USER root

# Update and upgrade packages
RUN apt-get update -qq && \
    apt-get -y upgrade && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install required packages
RUN apt-get update -qq && \
    mkdir -p /src && \
    apt-get -y --no-install-recommends install \
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
    jq && \
    apt-get -y --no-install-recommends install libbz2-dev liblzma-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

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
RUN apt-get update -qq && \
    apt-get -y install --no-install-recommends curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Drop back to the regular jenkins user - good practice
USER jenkins
