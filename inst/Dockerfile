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
