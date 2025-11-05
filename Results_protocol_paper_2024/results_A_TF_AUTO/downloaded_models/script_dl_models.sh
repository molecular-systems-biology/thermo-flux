#!/bin/bash

# Create a folder to store the models
download_folder="downloaded_models"
mkdir -p "$download_folder"

models=$(curl -s 'http://bigg.ucsd.edu/api/v2/models' | jq -r '.results[].bigg_id')

for model in $models
do
    echo "Downloading model: $model"
    curl -o "$download_folder/$model.sbml" "http://bigg.ucsd.edu/static/models/$model.xml"
done
