FROM rocker/verse
MAINTAINER Christian Panse <cp@fgcz.ethz.ch>
RUN install2.r --error \ 
  protViz \
  recmap\
  doParallel 
RUN install2.r --error \ 
  hexbin \
  rcdk 
RUN apt-get update \
  && apt-get install mono-complete vim less unzip r-base curl libxml2 -y
RUN echo "library(devtools); install_github('c-ruttkies/MetFragR/metfRag')" \
  | R --no-save 
RUN install2.r --error \ 
  lattice nlme
RUN echo "install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawDiag_0.0.33.tar.gz', repo=NULL)" \
  | R --no-save 
RUN apt-get install openjdk-8-jre:amd64 -y 
RUN install2.r --error testthat  
RUN R CMD javareconf
RUN echo "install.packages('rJava', type='source')" | R --no-save
RUN echo "install.packages('rcdk', type='source')" | R --no-save
RUN mkdir /tmp/UVPD
COPY . /tmp/UVPD/
# RUN cd /tmp && R CMD build UVPD && R CMD INSTALL /tmp/uvpd_*.tar.gz 
RUN cd /tmp && R CMD build UVPD --no-build-vignettes && R CMD INSTALL /tmp/uvpd_*.tar.gz 
RUN R --no-save < /tmp/UVPD/tests/testthat/test-install.R

# docker run -e PASSWORD=XXXX -v /srv/www/htdocs/p2722/:/srv/www/htdocs/p2722/  -p 8888:8787 1a5129ec532f 
