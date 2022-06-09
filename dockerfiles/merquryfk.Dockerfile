FROM ghcr.io/nbisweden/fastk:f18a4e6d2207539f7b84461daebc54530a9559b0 AS builder

ARG commit_id=8f3ab706e4cf4d7b7d1dfe5739859e3ebd26c494
RUN apt-get update && apt-get -y install \
    build-essential \
    git \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev

WORKDIR /opt
RUN git clone https://github.com/thegenemyers/MERQURY.FK.git merquryfk && \
    cd merquryfk && \
    git reset --hard ${commit_id} && \
    make all

FROM ghcr.io/nbisweden/fastk:f18a4e6d2207539f7b84461daebc54530a9559b0

ARG commit_id=8f3ab706e4cf4d7b7d1dfe5739859e3ebd26c494
LABEL description="FastK based version of Merqury " \
    container_author="Mahesh Binzer-Panchal" \
    version="${commit_id}" \
    author="Gene Myers" \
    repo="https://github.com/thegenemyers/MERQURY.FK"

COPY --from=builder /opt/merquryfk/ASMplot \
    /opt/merquryfk/CNplot \
    /opt/merquryfk/HAPmaker \
    /opt/merquryfk/HAPplot \
    /opt/merquryfk/KatComp \
    /opt/merquryfk/KatGC \
    /opt/merquryfk/MerquryFK \
    /opt/merquryfk/PloidyPlot \
    /usr/local/bin/

RUN apt-get update && apt-get -y install \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev

CMD [ "MerquryFK" ]