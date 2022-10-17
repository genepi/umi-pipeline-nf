FROM continuumio/miniconda
LABEL authors="Amstler Stephan" \
      email="amstler.stephan@i-med.ac.at" 

COPY environment.yml .
RUN \
   conda env update -n root -f environment.yml \
&& conda clean -a

# you only need to change <template> in this line:
# ENV PATH /opt/conda/envs/nf_umi_pipeline/bin:$PATH: