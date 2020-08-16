FROM rocker/shiny

MAINTAINER Pawel Zayakin "pawel@biomed.lu.lv"

RUN sed -i 's/main$/main contrib non-free/' /etc/apt/sources.list && \
  env DEBIAN_FRONTEND=noninteractive apt-get update && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends apt-utils && \
    env DEBIAN_FRONTEND=noninteractive apt-get upgrade -y --no-install-recommends && \
    env DEBIAN_FRONTEND=noninteractive apt-get dist-upgrade -y --no-install-recommends && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        git bowtie2 cutadapt samtools fastqc ncbi-blast+ python3-htseq && \
    apt-get autoremove -y && \
    apt-get autoclean -y
RUN R -e "chooseCRANmirror(graphics =FALSE,ind=1); \
          if (!requireNamespace('BiocManager')) install.packages('BiocManager'); \
          chooseBioCmirror(graphics =FALSE,ind=1); \
          BiocManager::install(c('gdata','shinydashboard','foreach','doMC','remotes'), ask=FALSE); \
          remotes::install_github('fbreitwieser/shinyFileTree', type = 'source')"
RUN mv /srv/shiny-server /srv/shiny-server.orig && \
    git clone https://github.com/zajakin/sRNAflow.git /srv/shiny-server && \
    mkdir /home/shiny/R && chmod 777 /home/shiny/R

# Estj: locales sudo wget gfortran libgfortran-10-dev libopenblas-dev \
#        libxml2-dev libjpeg-dev libcurl4-openssl-dev libssl-dev zlib1g-dev \
#        r-base-core git jupyter-notebook python3-pip \
#        bowtie bowtie2 cutadapt samtools fastqc ncbi-blast+ kraken2 python3-htseq rna-star \
#        fastp cnvkit seqtk picard-tools cufflinks bbmap trnascan-se trimmomatic radiant \
#        sortmerna bcftools gffread bedtools && \

#docker build  -t srnaflow shiny
#docker run -it --rm -p 3838:3838 -v `pwd`:/srv/shiny-server/www -v `pwd`/shinylog/:/var/log/shiny-server/ srnaflow
