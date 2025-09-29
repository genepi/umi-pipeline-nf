# Use micromamba minimal base image
FROM mambaorg/micromamba:1.5.8

LABEL authors="Amstler Stephan" \
      email="amstler.stephan@i-med.ac.at"

# Copy environment.yml into container
COPY environment.yml /tmp/environment.yml

# Install all dependencies into base environment, including wget
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba install -y -n base -c conda-forge wget && \
    micromamba clean --all --yes

# Ensure base environment binaries are on PATH (already default, but explicit)
ENV PATH=/opt/conda/bin:$PATH

# Make /opt writable by the default non-root user
USER root
RUN mkdir -p /opt && chown $MAMBA_USER:$MAMBA_USER /opt
USER $MAMBA_USER

# Set working directory
WORKDIR /opt

# Download mutserve jar as non-root
RUN wget https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc13.lpa/mutserve_LPA_adapted.jar
