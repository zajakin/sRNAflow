# sRNAflow


The goal of this project is to develop sRNAflow, a user-friendly bioinformatic tool with a web interface designed for the analysis of small RNAs obtained from biological fluids.
Pipeline includes filtering potential RNAs from reagents and environment, classifying small RNA types, managing small RNA annotation overlap, conducting differential expression assays, analysing isomiRs, approach to identify the sources of small RNAs within samples, alternative alignment-free analysis of RNA-seq data, featuring clustering and initial RNA source identification using BLAST.
Finalising results as a united report for all stages of workflow will be generated at the end of the complete process. Furthermore, the tool will be designed as user-friendly for inexperienced users.

To install sRNAflow on a server via Docker, execute the following command in your server terminal:
   docker pull ghcr.io/zajakin/srnaflow && chmod 777 . && \
   docker run -d -p 3838:3838 -v `pwd`:/srv/shiny-server/www ghcr.io/zajakin/srnaflow
   
After running, access the user interface in a web browser at HTTP://<your server name>:3838 (or another port if modified).
All uploads, databases, and analysis results are stored in the terminal’s current working folder. It is advised to create a local BLAST database using the command in the “Setup” tab.
Report files are accessible in the “Reports” tab. If you provide an email address and mail server, notifications and these report files can be sent for your convenience.
