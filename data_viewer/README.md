## ENCODE4 LR-RNA-seq data viewer
We've made it simple to visualize the ENCODE4 long-read RNA-seq data using both [Cerberus](https://github.com/mortazavilab/cerberus) and [Swan](https://github.com/mortazavilab/swan_vis) using this user-hosted data viewer. It runs using Docker so no software download is necessary.

#### Use instructions

0. Either clone this repository, or download the data download script.

Run one of the following in the bash terminal:

```bash
git clone git@github.com:fairliereese/paper_rnawg.git .
cd data_viewer
```

OR

```bash
wget https://raw.githubusercontent.com/fairliereese/paper_rnawg/master/data_viewer/download_data.sh
```

1. Download data and Docker image

This will automatically download the data files and put them into a `./data/` folder. Please do not move them. Download may take a while.

```bash
bash download_data.sh
```

2. Start the application

```bash
docker run --rm \
  -p 8501:8501 \
  -v $(pwd):/app \
  ghcr.io/fairliereese/encode4-viewer:latest
```

3. Open [http://localhost:8501](http://localhost:8501) in your browser. Happy gene hunting!
