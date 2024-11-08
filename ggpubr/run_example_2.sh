# example: generate a heatmap from csv
docker run \
    --rm \
    -v "$(pwd)/example_data":/mnt/data \
    rmillikin/ggpubr:latest \
    Rscript /mnt/data/plotHeatmap.R \
    --data /mnt/data/plotHeatmapData.csv \
    --aspect_ratio 1.0

# example: generate a volcano plot from csv
docker run \
    --rm \
    -v "$(pwd)/example_data":/mnt/data \
    rmillikin/ggpubr:latest \
    Rscript /mnt/data/plotVolcano.R \
    --data /mnt/data/plotVolcanoData.csv \
    --output /mnt/data/volcano.png \
    --pval_cutoff 0.05 \
    --l2fc_cutoff 1.0

# example: generate a confidence interval plot from csv

# TODO:
# violin plots
# bar plots
# rank (signal vs rank scatter)
# venn diagram
# could do scatter, line, box, etc.