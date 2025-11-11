# Global args
# MERQURY.FK: Release v1.1.3 + patch to save.cni file (09/2025)
# GENESCOPE.FK: Release v1.0 (full) (09/2023)
ARG merquryfk_commit_id=acef44f51ed5c431805682a42cc96616552b6cdb
ARG genescopefk_commit_id=380815c420f50171f9234a0fd1ff426b39829b91

FROM rocker/r-base:4.5.2 AS builder

ARG merquryfk_commit_id
ARG genescopefk_commit_id

RUN apt-get update && apt-get -y install \
    build-essential \
    git \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev

WORKDIR /opt

# Clone from patched version of MERQURY.FK
RUN git clone https://github.com/mahesh-panchal/MERQURY.FK.git merquryfk && \
    cd merquryfk && \
    git reset --hard ${merquryfk_commit_id} && \
    make all

FROM rocker/r-base:4.5.2

ARG merquryfk_commit_id
ARG genescopefk_commit_id

LABEL description="A container with MerquryFK and GeneScopeFK" \
    container_author="Mahesh Binzer-Panchal" \
    version="1.3" \
    merquryfk_commit_id="${merquryfk_commit_id}" \
    genescopefk_commit_id="${genescopefk_commit_id}" \
    tool_author="Gene Myers" \
    merquryfk_repo="https://github.com/thegenemyers/MERQURY.FK" \
    genescopefk_repo="https://github.com/thegenemyers/GENESCOPE.FK"

COPY --from=builder /opt/merquryfk/ASMplot \
    /opt/merquryfk/CNplot \
    /opt/merquryfk/HAPmaker \
    /opt/merquryfk/HAPplot \
    /opt/merquryfk/KatComp \
    /opt/merquryfk/KatGC \
    /opt/merquryfk/MerquryFK \
    /usr/local/bin/

RUN apt-get update && apt-get -y install \
    git \
    python3 \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    procps

WORKDIR /opt

RUN git clone https://github.com/thegenemyers/GENESCOPE.FK.git genescopefk && \
    cd genescopefk && \
    git reset --hard ${genescopefk_commit_id} && \
    Rscript -e 'install.packages(c("minpack.lm","argparse","ggplot2","scales","viridis","cowplot"))' && \
    Rscript -e 'install.packages(c("."), repos=NULL, type="source")' && \
    cp GeneScopeFK.R /usr/local/bin/

WORKDIR /

CMD [ "MerquryFK" ]
