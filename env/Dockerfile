FROM nfcore/base
LABEL authors="Adam Nunn" \
      email="adam.nunn@ecseq.com" \
      description="Docker image containing all requirements for the ecSeq/DNAseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# you only need to change <template> in this line:
ENV PATH /opt/conda/envs/<template>/bin:$PATH: