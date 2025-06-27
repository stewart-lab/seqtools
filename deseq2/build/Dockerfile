# Use a base image that contains R
FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libz-dev
#     imagemagick
#FROM stewartlab/deseq2:v1

# Install BiocManager and DESeq2
RUN R -e "install.packages(c('BiocManager', 'optparse'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('DESeq2', 'AnnotationDbi', 'org.Hs.eg.db', 'gprofiler2', 'GSVA', 'GSVABase', 'limma'))"

# Copy the data from ./src to /src in the image
RUN rm -rf ./src/
COPY ./src /src
RUN chmod +x src/*
# Set the entrypoint to run the R script using Rscript
# #["Rscript", "/src/deseq2.r"]
ENTRYPOINT ["Rscript", "/src/deseq2.r"]
#CMD ["R"]