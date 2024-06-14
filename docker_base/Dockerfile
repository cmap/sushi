FROM rocker/r-ver:4.0.3
MAINTAINER Anup Jonchhe <cmap-soft@broadinstitute.org>
LABEL base.mts.pipeline.clue.io.version="0.0.1"
LABEL base.mts.pipeline.clue.io.vendor="PRISM"

# Set environment variables for non-interactive installs
ENV DEBIAN_FRONTEND=noninteractive

# Update and install dependencies in a single RUN command to ensure proper order
RUN apt-get update -qq && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
      sudo \
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
      liblzma-dev && \
    # Ensure permissions for apt directories
    mkdir -p /var/cache/apt/archives/partial && \
    chown -R _apt:root /var/cache/apt/archives && \
    chmod -R 700 /var/cache/apt/archives/partial && \
    # Clean up apt cache to reduce image size
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

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
RUN apt-get update -qq && apt-get install -y --no-install-recommends curl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy and install additional sushi packages
COPY ./install_sushi_packages2.R /src/install_sushi_packages2.R
RUN Rscript /src/install_sushi_packages2.R

# Set the working directory
WORKDIR /usr/local/src/scripts

# Command to run your R script (replace 'your_script.R' with your main script)
# CMD ["Rscript", "your_script.R"]
