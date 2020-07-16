### Contact us

by E-mail [Brunner, Andrea](mailto: Andrea.Brunner@kwrwater.nl?SUBJECT=-> help uvpd shiny) or [cp@fgcz.ethz.ch](mailto:cp@fgcz.ethz.ch?SUBJECT=-> help uvpd shiny)

### Bug Reports

https://github.com/cpanse/uvpd/issues


### Run the  application

Install the R package

```{r}
pkgs <- c('shiny', 'ggplot2')
pkgs <- pkgs[(!pkgs %in% unique(installed.packages()[,'Package']))]
if(length(pkgs) > 0){install.packages(pkgs)}

install.packages('http://fgcz-ms.uzh.ch/~cpanse/UVPD/uvpd_0.0.9.tar.gz',repos=NULL)
```

run the shiny application from your computer

```{r}
shiny::runApp(file.path(system.file(package = 'uvpd'), 'shiny/stackedbarchart'))
```

### Useful bookmarks

- [Evaluation of food-relevant chemicals in the ToxCast high-throughput screening program](https://doi.org/10.1016/j.fct.2016.04.012)

- [Chemical Informatics Functionality in R](http://dx.doi.org/10.18637/jss.v018.i05)

- https://CRAN.R-project.org/package=rcdk

- ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt

- A new metaphor for projection-based visual analysis and data exploration Proc. SPIE
2007 | conference-paper [DOI: 10.1117/12.697879](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/6495/1/A-new-metaphor-for-projection-based-visual-analysis-and-data/10.1117/12.697879.short?SSO=1)

- https://CRAN.R-project.org/package=protViz

- [rawDiag - Brings Orbitrap Mass Spectrometry Data to Life; Multi-platform, Fast and Colorful R package. ](https://github.com/fgcz/rawDiag)
