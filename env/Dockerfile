FROM ubuntu:22.04

LABEL authors="Amstler Stephan" \
      email="amstler.stephan@i-med.ac.at"

COPY environment.yml .

RUN apt-get update && \
    apt-get install -y wget && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_23.9.0-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:${PATH}

RUN conda update -y conda && \
    conda env update -n root -f environment.yml && \
    conda clean --all

WORKDIR "/opt"
RUN wget https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc13.lpa/mutserve_LPA_adapted.jar
