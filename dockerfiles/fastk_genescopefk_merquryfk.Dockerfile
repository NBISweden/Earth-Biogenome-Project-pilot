FROM rocker/r-base:4.2.0 AS builder

ARG fastk_commit_id=f18a4e6d2207539f7b84461daebc54530a9559b0
# ARG merquryfk_commit_id=8f3ab706e4cf4d7b7d1dfe5739859e3ebd26c494
ARG merquryfk_commit_id=8ae344092df5dcaf83cfb7f90f662597a9b1fc61
ARG genescope_commit_id=380815c420f50171f9234a0fd1ff426b39829b91

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
# RUN git clone https://github.com/thegenemyers/MERQURY.FK.git merquryfk && \
RUN git clone https://github.com/mahesh-panchal/MERQURY.FK.git merquryfk && \
    cd merquryfk && \
    git reset --hard ${merquryfk_commit_id} && \
    make all

FROM rocker/r-base:4.2.0

ARG fastk_commit_id=f18a4e6d2207539f7b84461daebc54530a9559b0
ARG merquryfk_commit_id=8f3ab706e4cf4d7b7d1dfe5739859e3ebd26c494
ARG genescope_commit_id=380815c420f50171f9234a0fd1ff426b39829b91

LABEL description="A container with FastK, GeneScopeFK, and MerquryFK" \
    container_author="Mahesh Binzer-Panchal" \
    version="1.1" \
    fastk_commit_id="${fastk_commit_id}" \
    merquryfk_commit_id="${merquryfk_commit_id}" \
    genescope_commit_id="${genescope_commit_id}" \
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
    /opt/merquryfk/PloidyPlot \
    /usr/local/bin/

RUN apt-get update && apt-get -y install \
    git \
    python3 \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev

WORKDIR /opt
RUN git clone https://github.com/thegenemyers/GENESCOPE.FK.git genescopefk && \
    cd genescopefk && \
    git reset --hard ${genescope_commit_id} && \
    Rscript -e 'install.packages(c("minpack.lm","argparse","ggplot2","scales","viridis","cowplot"))' && \
    Rscript -e 'install.packages(c("."), repos=NULL, type="source")' && \
    cp GeneScopeFK.R /usr/local/bin/

WORKDIR /

CMD [ "FastK" ]