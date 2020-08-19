FROM rocker/shiny

MAINTAINER Pawel Zayakin "pawel@biomed.lu.lv"

RUN sed -i 's/main$/main contrib non-free/' /etc/apt/sources.list && \
  env DEBIAN_FRONTEND=noninteractive apt-get update && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends apt-utils && \
    env DEBIAN_FRONTEND=noninteractive apt-get upgrade -y --no-install-recommends && \
    env DEBIAN_FRONTEND=noninteractive apt-get dist-upgrade -y --no-install-recommends && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        git gawk unzip pigz \
        bowtie2 cutadapt samtools fastqc ncbi-blast+ ncbi-entrez-direct python3-htseq trnascan-se && \
    apt-get autoremove -y && \
    apt-get autoclean -y
RUN R -e "chooseCRANmirror(graphics =FALSE,ind=1); \
          if (!requireNamespace('BiocManager')) install.packages('BiocManager'); \
          chooseBioCmirror(graphics =FALSE,ind=1); \
          BiocManager::install(c('gdata','gplots','ggplot2','gridExtra','shinydashboard','DT','corrplot','edgeR','DESeq2', \
          'VennDiagram','GOstats','foreach','doMC','futile.logger','reticulate','shinyjs','sendmailR'), ask=FALSE)"
ARG DUMMY=unknown
RUN mv /srv/shiny-server /srv/shiny-server.orig && \
    DUMMY=${DUMMY} git clone https://github.com/zajakin/sRNAflow.git /srv/shiny-server && \
    git clone https://github.com/marbl/Krona.git /srv/shiny-server/Krona && /srv/shiny-server/Krona/KronaTools/install.pl && \
    git clone https://github.com/zajakin/ShortStack.git && \
    wget http://opengene.org/gencore/gencore -O /srv/shiny-server/bin/gencore && chmod a+x /srv/shiny-server/bin/gencore && \
    wget https://eda.polito.it/isomir-sea/isomiR-SEA_1.6_webpacket.zip -O tmp.zip && \
    unzip -oj tmp.zip isomiR-SEA_1.6_webpacket/isomiR-SEA_OS/Ubuntu_14_04_2LTS_x86_64/isomiR-SEA_1_6 && \
    mv isomiR-SEA_1_6 /srv/shiny-server/bin/isomiR-SEA && rm tmp.zip && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends bowtie libjson-perl
     # && mkdir /home/shiny/R && chmod 777 /home/shiny/R
# RUN R -e "chooseCRANmirror(graphics =FALSE,ind=1); chooseBioCmirror(graphics =FALSE,ind=1); \
#          BiocManager::install(c(), ask=FALSE)"
#                    remotes::install_github('fbreitwieser/shinyFileTree', type = 'source')"

# Estj: locales sudo wget gfortran libgfortran-10-dev libopenblas-dev \
#        libxml2-dev libjpeg-dev libcurl4-openssl-dev libssl-dev zlib1g-dev \
#        python3-pip kraken2 rna-star fastp cnvkit seqtk picard-tools trnascan-se sortmerna bcftools gffread bedtools radiant
#library(org.Hs.eg.db)

#docker build  -t srnaflow shiny
#chmod 777 .
#docker run -it --rm -p 3838:3838 -v `pwd`:/srv/shiny-server/www -v /tmp/shinylog/:/var/log/shiny-server/ srnaflow bash
