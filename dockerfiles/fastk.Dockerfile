FROM ubuntu:20.04 AS builder

ARG commit_id=f18a4e6d2207539f7b84461daebc54530a9559b0
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
    git reset --hard ${commit_id} && \
    make all

FROM ubuntu:20.04

ARG commit_id=f18a4e6d2207539f7b84461daebc54530a9559b0
LABEL description="A fast K-mer counter for high-fidelity shotgun datasets" \
    container_author="Mahesh Binzer-Panchal" \
    version="${commit_id}" \
    author="Gene Myers" \
    repo="https://github.com/thegenemyers/FASTK"

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

RUN apt-get update && apt-get -y install \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev

CMD [ "FastK" ]