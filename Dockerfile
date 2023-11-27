FROM openanalytics/r-base

MAINTAINER Jihun Kim "toujours209@gmail.com"

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libgsl-dev \
    libssh2-1-dev \
    libssl1.1 \
    libxml2-dev \
    build-essential \
    r-base-dev \
    pkg-config \
    cmake \
    libtiff5-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgdal-dev \
    && rm -rf /var/lib/apt/lists/*
    
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    unixodbc-dev \
    libmariadb-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('shiny', 'BiocManager', 'rmarkdown'), repos='https://cloud.r-project.org/')"
Run R -e "install.packages(c('seqinr', 'shinydashboard', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable'), repos = 'https://cloud.r-project.org/')"
RUN R -e "install.packages(c('DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'ape', 'compositions', 'stringr', 'caret', 'ggplot2', 'randomForest'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('data.table', 'xgboost', 'SHAPforxgboost', 'fontawesome', 'grid', 'ggplotify', 'remotes', 'doParallel', 'reshape2', 'fossil'), repos='https://cloud.r-project.org/')"
Run R -e "install.packages(c('proxy', 'ecodist', 'GUniFrac', 'picante', 'FSA', 'tibble', 'forestplot', 'VGAM', 'rgl', 'MiRKAT'), repos='https://cloud.r-project.org/')"

RUN R -e "BiocManager::install('phyloseq')"

RUN R -e "remotes::install_github('joey711/biomformat')"
RUN R -e "remotes::install_github('nik01010/dashboardthemes', force = TRUE)"
RUN R -e "remotes::install_github('zmjones/edarf', subdir = 'pkg')"
RUN R -e "remotes::install_github('Zhiwen-Owen-Jiang/MiRKATMC')"
RUN R -e "remotes::install_github('vegandevs/vegan')"

Run R -e "remotes::install_version('zCompositions', version = '1.4.0.1', repos = 'http://cran.us.r-project.org')"

RUN mkdir /root/app
COPY app /root/app
COPY Rprofile.site /usr/lib/R/etc/

COPY app/Data/biom.Rdata /root/app
COPY app/Data/otu.tab.txt /root/app
COPY app/Data/sam.dat.txt /root/app
COPY app/Data/tax.tab.txt /root/app
COPY app/Data/tree.tre /root/app

COPY app/www/MiMultiCat_Home_Img.png /root/app

COPY app/MiDataProc.Alpha.Diversity.R /root/app
COPY app/MiDataProc.Beta.Diversity.R /root/app
COPY app/MiDataProc.Data.Input.R /root/app
COPY app/MiDataProc.Data.Upload.R /root/app
COPY app/MiDataProc.ML.Models.R /root/app
COPY app/MiDataProc.ML.RF.R /root/app
COPY app/MiDataProc.ML.XGB.R /root/app
COPY app/MiDataProc.Taxa.R /root/app
COPY app/setSliderColor.R /root/app

EXPOSE 5512

CMD ["R", "-e", "shiny::runApp('/root/app')"]
