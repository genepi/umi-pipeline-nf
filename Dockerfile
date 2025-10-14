FROM continuumio/miniconda3:24.7.1-0

LABEL maintainer="Amstler Stephan <amstler.stephan@i-med.ac.at>" \
    version="v1.0.1"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    zlib1g-dev \
    libgomp1 \
    procps && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml

RUN conda install -y -n base -c conda-forge conda-libmamba-solver && \
    conda config --set solver libmamba && \
    conda env update -n base -f /tmp/environment.yml && \
    conda clean --all --yes

WORKDIR /opt

RUN wget -q https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc13.lpa/mutserve_LPA_adapted.jar

