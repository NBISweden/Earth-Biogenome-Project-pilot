# Global args
# FASTK: Release v1.1 + 12 unreleased commits (09/2025)
# MERQURY.FK: Release v1.1.3 + patch to save.cni file (09/2025)
# GENESCOPE.FK: Release v1.0 (full) (09/2023)
ARG fastk_commit_id=0e24fb45b71c4e14382ae1e1bc063bf66ea4e112
ARG merquryfk_commit_id=2ebdeb6a990f87edf2439acc353bba0cd34b07c7
ARG genescopefk_commit_id=380815c420f50171f9234a0fd1ff426b39829b91
ARG rbase_version=4.5.2

FROM rocker/r-base:${rbase_version} AS builder

ARG fastk_commit_id
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

RUN git clone https://github.com/thegenemyers/FASTK.git fastk && \
    cd fastk && \
    git reset --hard ${fastk_commit_id} && \
    make all

# Clone from patched version of MERQURY.FK
RUN git clone https://github.com/thegenemyers/MERQURY.FK.git merquryfk && \
    cd merquryfk && \
    git reset --hard ${merquryfk_commit_id} && \
    make all

FROM rocker/r-base:${rbase_version}

ARG fastk_commit_id
ARG merquryfk_commit_id
ARG genescopefk_commit_id

LABEL description="A container with Fastk, MerquryFK, and GeneScopeFK" \
    container_author="Mahesh Binzer-Panchal" \
    version="1.3" \
    fastk_commit_id="${fastk_commit_id}" \
    merquryfk_commit_id="${merquryfk_commit_id}" \
    genescopefk_commit_id="${genescopefk_commit_id}" \
    tool_author="Gene Myers" \
    fastk_repo="https://github.com/thegenemyers/FASTK" \
    merquryfk_repo="https://github.com/thegenemyers/MERQURY.FK" \
    genescopefk_repo="https://github.com/thegenemyers/GENESCOPE.FK"

COPY --from=builder /opt/fastk/FastK \
    /opt/fastk/Fastcp \
    /opt/fastk/Fastmv \
    /opt/fastk/Fastmerge \
    /opt/fastk/Fastrm \
    /opt/fastk/Haplex \
    /opt/fastk/Histex \
    /opt/fastk/Homex \
    /opt/fastk/Logex \
    /opt/fastk/Profex \
    /opt/fastk/Symmex \
    /opt/fastk/Tabex \
    /opt/fastk/Vennex \
    /usr/local/bin/

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

CMD [ "FastK" ]