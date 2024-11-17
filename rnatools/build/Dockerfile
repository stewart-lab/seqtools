FROM continuumio/miniconda3

# Install the required packages from the bioconda channel (fastp, rsem, fastqc)
# and automatically accept any prompts
RUN conda install -c bioconda fastqc fastp rsem star bioconductor-tximport -y

# Install the required R packages (BiocManager, optparse)
RUN R -e "install.packages(c('BiocManager', 'optparse'), repos='https://cloud.r-project.org')"

# Copy the data from ./src to /src in the image
COPY ./src /src

# Add entrypoint script
# COPY update_user.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /src/update_user.sh

# Set the entrypoint to the main.py script with unbuffered output
ENTRYPOINT ["python3", "-u", "/src/main.py"]