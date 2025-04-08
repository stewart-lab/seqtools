# create deseq2 input metadata (design, contrasts, counts)
# python create_deseq_metadata.py
# run --no-cache  when updating
# docker pull stewartlab/deseq2:v1
docker buildx build -t rmillikin/deseq2:latest .
# run deseq2
docker run \
    -i \
    -v /w5home/bmoore/sequencing/MSpurgeon/RNAseq/03.DEseq/:/mnt/data \
    rmillikin/deseq2:latest 
    # --design /mnt/data/design.csv \
    # --counts /mnt/data/counts.csv \
    # --contrasts /mnt/data/contrasts.csv \
    # --output_dir /mnt/data/