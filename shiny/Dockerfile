FROM rocker/shiny

MAINTAINER Pawel Zayakin "pawel@biomed.lu.lv"

RUN sed -i 's/main$/main contrib non-free/' /etc/apt/sources.list && \
  env DEBIAN_FRONTEND=noninteractive apt-get update && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends apt-utils && \
    env DEBIAN_FRONTEND=noninteractive apt-get upgrade -y --no-install-recommends && \
    env DEBIAN_FRONTEND=noninteractive apt-get dist-upgrade -y --no-install-recommends && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        git gawk unzip pigz libjson-perl curl python3-pip \
        bowtie bowtie2 cutadapt samtools fastqc ncbi-blast+ ncbi-entrez-direct python3-htseq trnascan-se && \
    apt-get autoremove -y && \
    apt-get autoclean -y
RUN R -e "chooseCRANmirror(graphics =FALSE,ind=1); \
          if (!requireNamespace('BiocManager')) install.packages('BiocManager'); \
          chooseBioCmirror(graphics =FALSE,ind=1); \
          BiocManager::install(c('gdata','gplots','ggplot2','gridExtra','shinydashboard','DT','corrplot','shinyjs','annotate'), ask=FALSE); \
          BiocManager::install(c('rtracklayer','VennDiagram','GOstats','foreach','doMC','futile.logger','sendmailR','openxlsx'), ask=FALSE); \
          BiocManager::install(c('XML','DESeq2'), ask=FALSE)"
# ARG DUMMY=unknown  # DUMMY=${DUMMY} git 
RUN mv /srv/shiny-server /srv/shiny-server.orig && \
    git clone https://github.com/zajakin/sRNAflow.git /srv/shiny-server && \
    git clone https://github.com/marbl/Krona.git /srv/shiny-server/Krona && cd /srv/shiny-server/Krona/KronaTools && ./install.pl && \
    git clone https://github.com/zajakin/ShortStack.git  /srv/shiny-server/ShortStack && \
    wget http://opengene.org/gencore/gencore -O /srv/shiny-server/bin/gencore && chmod a+x /srv/shiny-server/bin/gencore && \
    wget https://eda.polito.it/isomir-sea/isomiR-SEA_1.6_webpacket.zip -O tmp.zip && \
    unzip -oj tmp.zip isomiR-SEA_1.6_webpacket/isomiR-SEA_OS/Ubuntu_14_04_2LTS_x86_64/isomiR-SEA_1_6 && \
    mv isomiR-SEA_1_6 /srv/shiny-server/bin/isomiR-SEA && rm tmp.zip && \
    pip3 install multiqc mieaa
     # && mkdir /home/shiny/R && chmod 777 /home/shiny/R
# RUN env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends  && \
# RUN R -e "chooseCRANmirror(graphics =FALSE,ind=1); chooseBioCmirror(graphics =FALSE,ind=1); BiocManager::install(c('org.Hs.eg.db','edgeR','reticulate'), ask=FALSE)"
# remotes::install_github('fbreitwieser/shinyFileTree', type = 'source')"

#        libxml2-dev libjpeg-dev libcurl4-openssl-dev libssl-dev zlib1g-dev kraken2 rna-star fastp cnvkit seqtk picard-tools trnascan-se sortmerna bcftools gffread bedtools radiant

#chmod 777 . && docker pull zajakin/srnaflow && docker run -it --rm -p 3838:3838 -v `pwd`:/srv/shiny-server/www -v /tmp/shinylog/:/var/log/shiny-server/ zajakin/srnaflow
#docker build  -t srnaflow shiny && chmod 777 . && docker run -it --rm -p 3838:3838 -v `pwd`:/srv/shiny-server/www -v /tmp/shinylog/:/var/log/shiny-server/ srnaflow
