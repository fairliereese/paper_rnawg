FROM continuumio/miniconda3:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    ca-certificates \
    openssl \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# copy all the downloaded data files to the app
# COPY data/ /app/data

# conda stuff
# RUN conda config --set ssl_verify false # turn off ssl verification
# RUN conda config --remove channels defaults || true # remove defaults (repo.anaconda.com is blocked)
# RUN conda config --add channels conda-forge # add conda-forge as the only channel

# hard override conda config (kills defaults permanently)
RUN printf "\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
channel_priority: strict\n\
ssl_verify: false\n" > /opt/conda/.condarc

RUN conda config --set channel_priority strict
RUN conda install -n base -c conda-forge mamba

# create env
COPY environment.yml /app/environment.yml
RUN mamba create -n streamlit2 python=3.8.19
RUN conda run -n streamlit2 mamba env update -f /app/environment.yml

# pip stuff

###  github stuff

# first install git
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

# my versions of ternary, swan_vis, and cerberus
RUN conda run -n streamlit2 pip install git+http://github.com/fairliereese/python-ternary
RUN conda run -n streamlit2 pip install git+http://github.com/mortazavilab/swan_vis.git
RUN conda run -n streamlit2 pip install git+http://github.com/mortazavilab/cerberus.git


EXPOSE 8501

CMD ["conda", "run", "-n", "streamlit2", "streamlit", "run", "main.py", "--server.port=8501", "--server.address=0.0.0.0"]
