# dl data and put into data directory
mkdir data
wget https://zenodo.org/records/18299002/files/human_triplets.h5 -O data/human_triplets.h5
wget https://zenodo.org/records/18299002/files/mouse_triplets.h5 -O data/mouse_triplets.h5
wget https://zenodo.org/records/18299002/files/human_swan_graph.p -O data/human_swan.p
wget https://zenodo.org/records/18299002/files/mouse_swan_graph.p -O data/mouse_swan.p

# pull docker image
docker pull ghcr.io/fairliereese/encode4-viewer:latest
