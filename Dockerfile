FROM mambaorg/micromamba:1.5.8

LABEL authors="Amstler Stephan" \
      email="amstler.stephan@i-med.ac.at"

# Copy environment.yml
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

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