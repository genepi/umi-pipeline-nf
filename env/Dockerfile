FROM continuumio/miniconda
LABEL authors="Amstler Stephan" \
      email="amstler.stephan@i-med.ac.at" 

COPY environment.yml .
RUN \
   conda env update -n root -f environment.yml \
&& conda clean -a

#RUN mkdir /opt/mutserve
#WORKDIR "/opt/mutserve"
#RUN wget https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc15/mutserve.zip && \
#    unzip mutserve.zip
#ENV PATH="/opt/mutserve:${PATH}"

# you only need to change <template> in this line:
# ENV PATH /opt/conda/envs/nf_umi_pipeline/bin:$PATH:
