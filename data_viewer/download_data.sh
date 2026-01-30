# dl data and put into data director
mkdir data
wget https://zenodo.org/records/18415427/files/viewer_data.tar.gz -O data/viewer_data.tar.gz
tar -xzf data/viewer_data.tar.gz -C data/

# pull docker image
docker pull ghcr.io/fairliereese/encode4-viewer:latest
