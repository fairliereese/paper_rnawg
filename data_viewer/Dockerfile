# ---------- build stage ----------
FROM continuumio/miniconda3:latest AS builder

# --------------------------
# System dependencies for building packages
# --------------------------
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    ca-certificates \
    openssl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /wheels

# --------------------------
# Conda configuration (SSL off, channels)
# --------------------------
RUN printf "\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
channel_priority: strict\n\
ssl_verify: false\n" > /opt/conda/.condarc

ENV CONDA_SSL_VERIFY=no

# --------------------------
# Create Python environment & install PyTables/HDF5 via conda
# --------------------------
RUN conda create -y -n app_env python=3.8
RUN conda run -n app_env conda install -y -c conda-forge pytables=3.8.0 hdf5

# Upgrade pip & wheel inside the conda env
RUN conda run -n app_env pip install --upgrade pip wheel

# --------------------------
# Build wheels for other Python packages
# --------------------------
RUN conda run -n app_env pip wheel \
    --wheel-dir=/wheels \
    streamlit==1.40.1 \
    numpy==1.24.4 \
    pandas==1.5.3 \
    anndata==0.9.2 \
    scanpy==1.9.8 \
    scipy==1.10.1 \
    pyranges==0.1.4 \
    h5py==3.11.0 \
    plotly \
    git+https://github.com/fairliereese/python-ternary \
    git+https://github.com/mortazavilab/swan_vis.git \
    git+https://github.com/mortazavilab/cerberus.git

# ------------ runtime stage---------------
FROM continuumio/miniconda3:latest

# Copy wheels from builder stage
COPY --from=builder /wheels /wheels


RUN apt-get update && apt-get install -y \
    ca-certificates \
    openssl \
    libblas-dev \
    liblapack-dev \
    gfortran \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy full conda environment (with PyTables, pip packages, etc.)
COPY --from=builder /opt/conda/envs/app_env /opt/conda/envs/app_env

# Conda SSL config for runtime
RUN printf "\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
channel_priority: strict\n\
ssl_verify: false\n" > /opt/conda/.condarc

ENV CONDA_SSL_VERIFY=no

# Install wheels into app_env
RUN conda run -n app_env pip install --no-cache-dir /wheels/*


# Activate environment in CMD
SHELL ["conda", "run", "-n", "app_env", "/bin/bash", "-c"]

EXPOSE 8501
CMD ["streamlit", "run", "main.py", "--server.port=8501", "--server.address=0.0.0.0"]
