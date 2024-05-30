FROM docker.io/library/r-base:4.4.0

WORKDIR /workspace

COPY install_packages.R /workspace

RUN apt-get update \
    && apt-get install -yq libcurl4-openssl-dev libfontconfig1-dev libfreetype6-dev \
    libfribidi-dev libharfbuzz-dev libhdf5-dev libjpeg-dev libpng-dev libssl-dev \
    libtiff5-dev libxml2-dev \
    && Rscript install_packages.R \
    && rm -rf /tmp/* \
    && rm -rf /var/cache/apt/* \
    && rm -rf /var/tmp/*

