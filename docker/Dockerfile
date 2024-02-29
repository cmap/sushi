FROM prismcmap/base-clue-pseq
#FROM prismcmap/base-clue-pseq-test
MAINTAINER Anup Jonchhe <cmap-soft@broadinstitute.org>
LABEL clue.sushi.pipeline.clue.io.version="0.0.1"
LABEL clue.sushi.pipeline.clue.io.vendor="PRISM/CMAP"

RUN apt-get -y update

RUN apt-get install -y curl \
    bzip2 \
    git

RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/miniconda/bin:${PATH}

RUN conda create -y --name prism -c bioconda python=3.6 numpy pandas yaml h5py requests setuptools argparse pathlib pyyaml jinja2 simplejson

RUN mkdir -p /cmap/bin && \
mkdir -p /cmap/bin/sushi/

#COPY ../R /cmap/bin/sushi/R
#COPY ../scripts /cmap/bin/sushi/scripts
COPY ./ /cmap/bin/sushi/
COPY ./docker/install_packages.R /src/install_packages.R

RUN Rscript /src/install_packages.R
RUN Rscript -e "setwd('/cmap/bin/sushi/');getwd();devtools::install()"
COPY ./docker/aws_batch.sh /clue/bin/aws_batch

WORKDIR /
ENV PATH /clue/bin:$PATH
RUN ["chmod","-R", "+x", "/cmap/bin"]
#RUN ["chmod","-R", "+x", "/cmap/sushi"]
ENTRYPOINT ["aws_batch"]

CMD ["-help"]
