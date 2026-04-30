# create deseq2 input metadata (design, contrasts, counts)
# python create_deseq_metadata.py
# run --no-cache  when updating
# docker pull stewartlab/deseq2:v1 #rmillikin/deseq2:latest
#docker buildx build -t rmillikin/deseq2:latest .
# run deseq2
docker run \
    -ti --rm \
    -v /w5home/bmoore/TTN_HFpEF_risk/:/mnt/data \
    rmillikin/deseq2:latest \
    --design /mnt/data/design.csv \
    --counts /mnt/data/HFpEF_study_BAMs.filtered_reassigned_counts_high-lowrisk_techmeans.csv \
    --contrasts /mnt/data/contrasts.csv \
    --output_dir /mnt/data/