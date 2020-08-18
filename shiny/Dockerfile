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
          BiocManager::install(c('gdata','gplots','ggplot2','gridExtra','shinydashboard','DT','corrplot','edgeR','DESeq2','VennDiagram','GOstats','foreach','doMC','futile.logger','reticulate','shinyjs'), ask=FALSE)"
ARG DUMMY=unknown
RUN mv /srv/shiny-server /srv/shiny-server.orig && \
    DUMMY=${DUMMY} git clone https://github.com/zajakin/sRNAflow.git /srv/shiny-server && \
    git clone https://github.com/marbl/Krona.git /Krona && /Krona/install.pl && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends gawk radiant trnascan-se bowtie
     # && mkdir /home/shiny/R && chmod 777 /home/shiny/R
RUN R -e "chooseCRANmirror(graphics =FALSE,ind=1); chooseBioCmirror(graphics =FALSE,ind=1); \
         BiocManager::install(c('sendmailR'), ask=FALSE)"
#                    remotes::install_github('fbreitwieser/shinyFileTree', type = 'source')"

# Estj: locales sudo wget gfortran libgfortran-10-dev libopenblas-dev \
#        libxml2-dev libjpeg-dev libcurl4-openssl-dev libssl-dev zlib1g-dev \
#        python3-pip kraken2 rna-star fastp cnvkit seqtk picard-tools trnascan-se sortmerna bcftools gffread bedtools
#library(org.Hs.eg.db)

#docker build  -t srnaflow shiny
#chmod 777 .
#docker run -it --rm -p 3838:3838 -v `pwd`:/srv/shiny-server/www -v /tmp/shinylog/:/var/log/shiny-server/ srnaflow bash
