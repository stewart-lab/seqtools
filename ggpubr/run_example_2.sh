# an example where you can run an R script in the container for reproducible plots.
docker run --rm -v "$(pwd)/example_data":/mnt/data rmillikin/ggpubr:latest Rscript /mnt/data/exampleScript.R