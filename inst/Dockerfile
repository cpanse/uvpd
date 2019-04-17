FROM rocker/verse
MAINTAINER Christian Panse <cp@fgcz.ethz.ch>
RUN install2.r --error \ 
  protViz \
  recmap\
  doParallel 
RUN install2.r --error \ 
  hexbin \
  rcdk 
RUN echo "install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawDiag_0.0.33.tar.gz', repo=NULL)" \
  | R --no-save 
RUN echo "library(devtools); install_github('c-ruttkies/MetFragR/metfRag')" \
  | R --no-save 
RUN mkdir /tmp/UVPD
COPY . /tmp/UVPD/
RUN cd /tmp && find UVPD && R CMD build UVPD --no-build-vignettes && R CMD INSTALL /tmp/uvpd_*.tar.gz
