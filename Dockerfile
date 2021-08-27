# docker build --no-cache -t singjust/dialignr:2.0.0 .
# docker push singjust/dialignr:2.0.0

FROM rocker/r-ver:4.1.0

#######################
##      System
####################### 

RUN apt-get update && \
    apt-get install -y zlib1g-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libssl-dev && \
    apt-get install -y libnetcdf-dev && \
    apt-get install -y libxml2-dev && \
    apt-get install -y libglpk-dev

# RUN R -e "if(!requireNamespace('BiocManager', quietly = TRUE)){ install.packages('BiocManager')}; BiocManager::install('DIAlignR')"

RUN R -e "if(!requireNamespace('devtools', quietly = TRUE)){ install.packages('devtools')}; devtools::install_github('singjc/DIAlignR', ref='feature/docker')"

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)){ install.packages('BiocManager')}; BiocManager::install('BiocParallel')"

COPY Rscript/alignTargetedRuns_cli.R /alignTargetedRuns_cli.R

ENTRYPOINT ["Rscript", "alignTargetedRuns_cli.R"]
CMD ["--help"]
