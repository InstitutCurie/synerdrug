## Emacs, make this -*- mode: sh; -*-

FROM rocker/shiny:latest

ENV http_proxy="http://www-cache.curie.fr:3128"
ENV https_proxy="http://www-cache.curie.fr:3128"

RUN R -e "install.packages(c('shinydashboard', 'shinyBS', 'ggplot2', 'DT', 'RColorBrewer', 'reshape2', 'dplyr', 'mgcv', 'magrittr', 'minpack.lm'), repos='http://cran.rstudio.com/')"


CMD ["R -e 'require(synerdrug); runSynerDrug()'"]
