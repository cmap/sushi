FROM prismcmap/base-clue-pseq
#FROM prismcmap/base-clue-pseq-test
MAINTAINER Anup Jonchhe <cmap-soft@broadinstitute.org>
LABEL clue.sushi.pipeline.clue.io.version="0.0.1"
LABEL clue.sushi.pipeline.clue.io.vendor="PRISM/CMAP"

RUN apt-get -y update
RUN apt-get -y install git

RUN mkdir -p /cmap/bin && \
mkdir -p /cmap/bin/sushi/

RUN cd /cmap/bin/ && \
git clone https://github.com/cmap/sushi.git && \
cd /

COPY ./aws_batch.sh /clue/bin/aws_batch

WORKDIR /
ENV PATH /clue/bin:$PATH
RUN ["chmod","-R", "+x", "/cmap/bin"]
#RUN ["chmod","-R", "+x", "/cmap/sushi"]
ENTRYPOINT ["aws_batch"]

CMD ["-help"]