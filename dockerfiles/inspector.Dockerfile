# https://hub.docker.com/r/mambaorg/micromamba
FROM mambaorg/micromamba:0.15.2

SHELL ["/bin/bash", "-c"]

# Add metadata to the container using labels.
LABEL description="Temporary Inspector container" \
    author="Mahesh Binzer-Panchal" \
    version="1.0.0"

# Change user to root to install things (set to micromamba in base image)
USER root

# Update and install software dependencies
# with APT (Advanced Packaging  Tool)
RUN apt-get update && \
    apt-get install -y procps ghostscript git

# Install tools using the mamba package manager
RUN micromamba install -y -n base \
        conda-forge::python=2.7.15 \
        bioconda::pysam=0.16.0.1 \
        bioconda::minimap2=2.15 \
        bioconda::samtools=1.9 \
        bioconda::flye=2.8.3 \
        anaconda::statsmodels=0.10.1 && \
    micromamba clean --all --yes

# Install Inspector in opt
WORKDIR /opt
RUN git clone --depth 1 https://github.com/Maggi-Chen/Inspector.git

ENV PATH="/opt/Inspector:${PATH}"

CMD [ "inspector.py" ]
