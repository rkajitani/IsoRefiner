FROM continuumio/miniconda3
RUN conda install -y -c conda-forge -c bioconda isorefiner
