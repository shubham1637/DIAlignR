# docker build --no-cache -t singjust/dialignr:2.0.0 .
# docker push singjust/dialignr:2.0.0

FROM rocker/r-ver:3.4.4

RUN R -e "if(!requireNamespace('BiocManager', quietly = TRUE)){ install.packages('BiocManager')}; BiocManager::install('DIAlignR')"

