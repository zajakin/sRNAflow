FROM rocker/shiny

MAINTAINER Pawel Zayakin "pawel@biomed.lu.lv"

RUN sed -i 's/main$/main contrib non-free/' /etc/apt/sources.list && \
  env DEBIAN_FRONTEND=noninteractive apt-get update && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends apt-utils && \
    env DEBIAN_FRONTEND=noninteractive apt-get upgrade -y --no-install-recommends && \
    env DEBIAN_FRONTEND=noninteractive apt-get dist-upgrade -y --no-install-recommends && \
    env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        git gawk unzip pigz libjson-perl curl python3-pip libxml2-dev zlib1g-dev libbz2-dev liblzma-dev libssl-dev \
        bowtie bowtie2 cutadapt samtools fastqc ncbi-blast+ ncbi-entrez-direct python3-htseq trnascan-se && \
    apt-get autoremove -y && \
    apt-get autoclean -y
RUN R -e "chooseCRANmirror(graphics =FALSE,ind=1); \
          if (!requireNamespace('BiocManager')) install.packages('BiocManager'); \
          chooseBioCmirror(graphics =FALSE,ind=1); \
          BiocManager::install(c('gdata','gplots','ggplot2','gridExtra','shinydashboard','DT','corrplot','shinyjs','annotate'), ask=FALSE); \
          BiocManager::install(c('rtracklayer','VennDiagram','GOstats','foreach','doMC','futile.logger','sendmailR','openxlsx'), ask=FALSE); \
          BiocManager::install(c('XML','DESeq2'), ask=FALSE)"
# ARG DUMMY=unknown  
# RUN DUMMY=${DUMMY} ls
RUN mv /srv/shiny-server /srv/shiny-server.orig && mkdir -p /srv/shiny-server/bin && \
    git clone https://github.com/marbl/Krona.git /srv/shiny-server/Krona && cd /srv/shiny-server/Krona/KronaTools && ./install.pl && \
    git clone https://github.com/zajakin/ShortStack.git  /srv/shiny-server/ShortStack && \
    wget http://opengene.org/gencore/gencore -O /srv/shiny-server/bin/gencore && chmod a+x /srv/shiny-server/bin/gencore && \
    wget https://eda.polito.it/isomir-sea/isomiR-SEA_1.6_webpacket.zip -O tmp.zip && \
    unzip -oj tmp.zip isomiR-SEA_1.6_webpacket/isomiR-SEA_OS/Ubuntu_14_04_2LTS_x86_64/isomiR-SEA_1_6 && \
    mv isomiR-SEA_1_6 /srv/shiny-server/bin/isomiR-SEA && rm tmp.zip && \
    pip3 install multiqc mieaa
RUN pip3 install tensorflow keras && git clone https://github.com/tjgu/miTAR.git /srv/shiny-server/miTAR && \
    sed -i 's/tf.set_random_seed(sdnum)/tf.random.set_seed(sdnum)/p' /srv/shiny-server/miTAR/predict_multimiRmultimRNA.py
COPY app.R /srv/shiny-server/app.R
COPY LICENSE /srv/shiny-server/LICENSE
COPY README.md /srv/shiny-server/README.md
COPY bin /srv/shiny-server/bin/
COPY shiny /srv/shiny-server/shiny/
COPY gtf_biotypes /srv/shiny-server/gtf_biotypes/
    # git clone https://github.com/zajakin/sRNAflow.git /srv/shiny-server && \
     # && mkdir -m 777 /home/shiny/R
# RUN wget http://cbio.mskcc.org/miRNA2003/src1.9/binaries/miRanda-1.9-i686-linux-gnu.tar.gz
# RUN wget https://github.com/cbiagii/target-prediction/blob/master/miRanda-aug2010.tar.gz?raw=true -O miRanda-aug2010.tar.gz && tar zxvf miRanda-aug2010.tar.gz && cd miRanda-3.3a && ./configure && make
# RUN env DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends  && \
# RUN R -e "chooseCRANmirror(graphics =FALSE,ind=1); chooseBioCmirror(graphics =FALSE,ind=1); BiocManager::install(c('org.Hs.eg.db','edgeR','reticulate'), ask=FALSE)"
# remotes::install_github('fbreitwieser/shinyFileTree', type = 'source')"

#        libjpeg-dev libcurl4-openssl-dev kraken2 rna-star fastp cnvkit seqtk picard-tools sortmerna bcftools gffread bedtools radiant

#chmod 777 . && docker pull zajakin/srnaflow && docker run -it --rm -p 3838:3838 -v `pwd`:/srv/shiny-server/www -v /tmp/shinylog/:/var/log/shiny-server/ zajakin/srnaflow
#docker build  -t srnaflow . && chmod 777 . && docker run -it --rm -p 3838:3838 -v `pwd`:/srv/shiny-server/www -v /tmp/shinylog/:/var/log/shiny-server/ srnaflow
