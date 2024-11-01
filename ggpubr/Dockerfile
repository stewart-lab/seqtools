# Use a base image that contains R, RStudio
FROM rocker/verse:4.4.1

RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libz-dev

# Install ggpubr and optparse
RUN R -e "install.packages(c('ggpubr', 'optparse'), repos='https://cloud.r-project.org')"