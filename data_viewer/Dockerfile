

# ---------- build stage ----------
FROM python:3.8-slim AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /wheels

RUN pip install --upgrade pip wheel

# install everything into wheels
RUN pip wheel \
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


# ---------- runtime stage ----------
FROM python:3.8-slim

RUN apt-get update && apt-get install -y \
    ca-certificates \
    openssl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /wheels /wheels

RUN pip install --no-cache-dir /wheels/*

EXPOSE 8501
CMD ["streamlit", "run", "main.py", "--server.port=8501", "--server.address=0.0.0.0"]
