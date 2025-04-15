#docker pull rmillikin/rnatools:latest
#docker pull rmillikin/deseq2:latest

# align to human genome
#     -d \
#   --rm \
docker build -t rmillikin/rnatools:latest ./ 
docker run \
    -d \
    --rm \
    -v /w5home/bmoore/sequencing:/w5home/bmoore/sequencing \
    rmillikin/rnatools:latest \
    --fastq_dir /w5home/bmoore/sequencing/MSpurgeon/RNAseq/04.Reruns/01.RawData \
    --reference_genome_dir /w5home/bmoore/sequencing/Homo_sapiens_CRCh38 \
    --output_dir /w5home/bmoore/sequencing/MSpurgeon/RNAseq/04.Reruns/02.alignment \
    --resume

# align to mouse genome
# docker run \
#     -d \
#     --rm \
#     -v /w5home/rmillikin/sequencing:/w5home/rmillikin/sequencing \
#     rmillikin/rnatools:latest \
#     --fastq_dir /w5home/rmillikin/sequencing/kratz/2024-10-16_Pan33_Pan64_RNASeq/analysis/000_fastq \
#     --reference_genome_dir /w5home/rmillikin/sequencing/reference_genomes/Mus_musculus_GRCm39 \
#     --output_dir /w5home/rmillikin/sequencing/kratz/2024-10-16_Pan33_Pan64_RNASeq/analysis \
#     --resume

# create deseq2 input metadata (design, contrasts, counts)
# python create_deseq_metadata.py

# run deseq2
# docker run \
#     --rm -v "$(pwd)/analysis/2024-11-05_15-23-42___rnatools":/mnt/data \
#     rmillikin/deseq2:latest \
#     --design /mnt/data/006_deseq_input/design.csv \
#     --counts /mnt/data/006_deseq_input/counts.csv \
#     --contrasts /mnt/data/006_deseq_input/contrasts.csv \
#     --output_dir /mnt/data/007_deseq_output