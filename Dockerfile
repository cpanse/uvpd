FROM rocker/verse
MAINTAINER Christian Panse <cp@fgcz.ethz.ch>
LABEL version="uvpd 0.0.1"
LABEL description="docker image contains repro code snippets for Project 2722 - Benchmark on Ultra HRMS in combination with UVPD fragmentation for enhanced structural identification of organic micropollutants"
RUN apt-get update \
  && apt-get install mono-complete vim less unzip r-base curl libxml2 openjdk-8-jre:amd64 -y
RUN echo "install.packages('rJava', type='source')" | R --no-save \
  && R CMD javareconf
RUN install2.r --error \ 
  lattice \
  protViz \
  doParallel \
  hexbin \
  rcdk \
  testthat \
  deisotoper \
  nlme
RUN echo "library(devtools); install_github('c-ruttkies/MetFragR/metfRag')" \
  | R --no-save 
RUN echo "install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawDiag_0.0.33.tar.gz', repo=NULL)" \
  | R --no-save 
#RUN echo "install.packages('rcdk', type='source')" | R --no-save
RUN mkdir /tmp/UVPD
COPY . /tmp/UVPD/
# RUN cd /tmp && R CMD build UVPD && R CMD INSTALL /tmp/uvpd_*.tar.gz 
RUN cd /tmp && R CMD build UVPD --no-build-vignettes && R CMD INSTALL /tmp/uvpd_*.tar.gz 
RUN R --no-save < /tmp/UVPD/tests/testthat/test-install.R
EXPOSE 8888:8787 


# docker run --name p2772_uvpd -e PASSWORD=XXXX -v /srv/www/htdocs/p2722/:/srv/www/htdocs/p2722/ -p 8888:8787 d3dc766214da
