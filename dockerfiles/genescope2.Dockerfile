from rocker/r-base:4.2.0

ARG commit_id=380815c420f50171f9234a0fd1ff426b39829b91
LABEL description="A derivative of GenomeScope2.0 modified to work with FastK" \
    container_author="Mahesh Binzer-Panchal" \
    version="${commit_id}" \
    author="Gene Myers" \
    repo="https://github.com/thegenemyers/GENESCOPE.FK"

RUN apt-get update && apt-get -y install \
    git=1:2.35.1-1 \
    python3=3.10.4-1+b1

WORKDIR /opt
RUN git clone https://github.com/thegenemyers/GENESCOPE.FK.git genescopefk && \
    cd genescopefk && \
    git reset --hard ${commit_id} && \
    Rscript -e 'install.packages(c("minpack.lm","argparse"))' && \
    Rscript -e 'install.packages(c("."), repos=NULL, type="source")' && \
    cp GeneScopeFK.R /usr/local/bin/

CMD [ "GeneScopeFK.R" ]