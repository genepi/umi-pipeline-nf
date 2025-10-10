FROM mambaorg/micromamba:1.5.8

LABEL authors="Amstler Stephan" \
      email="amstler.stephan@i-med.ac.at"

# Copy environment.yml
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Install pip-only dependencies in a separate layer
RUN micromamba run -n base pip install \
    matplotlib==3.7.2 \
    networkx==3.1 \
    catfishq==1.4.0 \
    numpy==1.23.4 \
    pysam==0.19.1 \
    pyfastx==1.1.0 \
    pandas==1.5.0 \
    edlib==1.3.9

# Installing software using system package manager
USER root

# Install software
RUN apt-get update && \
    apt-get install -y \
    wget \
    zlib1g-dev \
    libgomp1 \
    procps && \
    rm -rf /var/lib/apt/lists/*

# Working directory
WORKDIR /opt

# Download mutserve jar
RUN wget https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc13.lpa/mutserve_LPA_adapted.jar
