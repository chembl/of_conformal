# OpenFaaS Conformal Prediction models

Re-implementation of https://doi.org/10.1186/s13321-018-0325-4 with [LightGBM](https://lightgbm.readthedocs.io/en/latest/).


Predicting only models with CCR ((sensitivity + specificity) / 2) >= 0.85

# Installation

## Build the image:
```
faas-cli build -f of_conformal.yml
```

## Deploy it to OpenFaaS
```
faas-cli deploy -f of_conformal.yml
```

## Run it locally (image available in [DockerHub](https://hub.docker.com/r/chembl/mcp))

Since it needs to load the models into memory it may take few minutes to start properly replying.

```
docker run -p 8080:8080 chembl/mcp

curl -X POST -H 'Accept: */*' -H 'Content-Type: application/json' -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}' http://127.0.0.1:8080
```
