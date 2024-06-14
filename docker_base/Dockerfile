FROM rocker/r-ver:4.0.3
MAINTAINER Anup Jonchhe <cmap-soft@broadinstitute.org>
LABEL base.mts.pipeline.clue.io.version="0.0.1"
LABEL base.mts.pipeline.clue.io.vendor="PRISM"

USER root

# Update and install dependencies
RUN apt-get update -qq && \
    apt-get -y upgrade && \
    apt-get install -y sudo && \
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
      jq \
      libbz2-dev \
      liblzma-dev

# Create a custom library directory with appropriate permissions
RUN mkdir -p /usr/local/lib/R/custom-library && \
    chmod 777 /usr/local/lib/R/custom-library

# Set the R_LIBS_USER environment variable to the custom library path
ENV R_LIBS_USER=/usr/local/lib/R/custom-library

# Set TMPDIR to a writable directory
ENV TMPDIR=/tmp

# Copy and install R packages
COPY ./install_packages.R /src/install_packages.R
RUN Rscript /src/install_packages.R

# Copy and install sushi packages
COPY ./install_sushi_packages.R /src/install_sushi_packages.R
RUN Rscript /src/install_sushi_packages.R

# Install curl
RUN apt-get -y update && apt-get -y install curl

# Copy and install additional sushi packages
COPY ./install_sushi_packages2.R /src/install_sushi_packages2.R
RUN Rscript /src/install_sushi_packages2.R

# Set the working directory
WORKDIR /usr/local/src/scripts

# Command to run your R script (replace 'your_script.R' with your main script)
# CMD ["Rscript", "your_script.R"]
