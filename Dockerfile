# syntax=docker/dockerfile:1

FROM mambaorg/micromamba:1.5.8

LABEL authors="Amstler Stephan" \
      email="amstler.stephan@i-med.ac.at"

# Install everything via conda-forge/bioconda, except catfishq (pip)
RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    python>=3.8 \
    medaka=2.0.1 \
    minimap2=2.24 \
    samtools=1.15.1 \
    seqtk=1.3 \
    lofreq=2.1.5 \
    freebayes=1.3.2 \
    vcflib=1.0.0 \
    bedtools=2.30.0 \
    vsearch=2.21.2 \
    openjdk=11.0.9 \
    unzip=6.0 \
    matplotlib=3.7.2 \
    networkx=3.1 \
    numpy=1.23.4 \
    pysam=0.19.1 \
    pyfastx=1.1.0 \
    pandas=1.5.0 \
    edlib=1.3.9 \
    pip && \
    micromamba clean --all --yes

# Only pip-install catfishq (not on conda)
RUN micromamba run -n base pip install catfishq==1.4.0

# Install Linux tools
USER root
RUN apt-get update && \
    apt-get install -y \
    wget \
    zlib1g-dev \
    libgomp1 \
    procps && \
    rm -rf /var/lib/apt/lists/*

# Working directory
WORKDIR /opt

# Download mutserve JAR
RUN wget https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc13.lpa/mutserve_LPA_adapted.jar
