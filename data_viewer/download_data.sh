# dl data and put into data directory
mkdir data
wget link https://zenodo.org/records/18299002/files/human_triplets.h5 -o data/human_triplets.h5
wget link https://zenodo.org/records/18299002/files/mouse_triplets.h5 -o data/mouse_triplets.h5
wget link https://zenodo.org/records/18299002/files/human_swan_graph.p -o data/human_swan.p
wget link https://zenodo.org/records/18299002/files/mouse_swan_graph.p -o data/mouse_swan.p

# pull docker image
docker pull ghcr.io/fairliereese/encode4-viewer:latest
