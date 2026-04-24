FROM rocker/shiny

RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    build-essential \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libglpk40 \
    libproj-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN R -e 'install.packages(c("DT","data.table","dplyr","tidyr", \
    "pheatmap","ggplot2","cowplot","corrplot","plotly","stringr", \
    "BiocManager","devtools", "enrichR"), dependencies=TRUE)'

# Install Bioconductor core deps first
RUN R -e 'BiocManager::install(c("DESeq2","edgeR","AnnotationDbi","org.Hs.eg.db"), ask=FALSE)'

# Install PCAtools (after dependencies)
RUN R -e 'BiocManager::install("PCAtools", ask=FALSE)'

RUN rm -rf /srv/shiny-server/*
WORKDIR /srv/shiny-server/
COPY ./ ./