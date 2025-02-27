# Use a base image that contains R
FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libz-dev

# Install BiocManager and DESeq2
RUN R -e "install.packages(c('BiocManager', 'optparse'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('DESeq2', 'AnnotationDbi', 'org.Hs.eg.db', 'gprofiler2', 'GSVA', 'GSVABase', 'limma'))"

# Copy the data from ./src to /src in the image
COPY ./src /src

# Set the entrypoint to run the R script using Rscript
ENTRYPOINT ["Rscript", "/src/deseq2.r"]
