# Bayesian-evolvability-NCII

Bayesian analysis of *Molgula occidentalis* (ascidian) larval traits from a North Carolina II (NCII) breeding design. <br />

The repository contains R code and data to produce the analyses and plots in the manuscipt Powell JA, Archer FI, Burgess SC. Evolvability of larval dispersal traits in a marine invertebrate. <br />

### R Package Dependencies
`rjags` and local installation of JAGS (https://mcmc-jags.sourceforge.io/) <br />
`evolvability` <br />
`QGglmm` <br />
`swfscMisc` <br />
`abind` <br />
`tidyverse` <br />

### Files

`Moccidentalis All Traits.csv` raw data file. <br />
`head_tail_data.rds` formatted data for Model I and II in RDS format (see `?readRDS`).  <br />
`hatch_settle.rds` formatted data for Models III and IV in RDS format (see `?readRDS`).  <br />
`Model_##.R` scripts for running Models I - IV.  <br />
`Model_##.Rmd` Rmarkdown file for summarizing posterior output from model runs.  <br />
`Figure_Model ##.R` scripts for making figures in manuscript.  <br />
`0_misc_funcs.R` miscellaneous functions used in summary and figure making files.  <br />
`1_load_and_format.R` load raw data and create formatted RDS files for use in models.  <br />

### Contact
For questions about the code, please contact [Eric Archer](https://github.com/EricArcher). <br />
For bug reports, please leave an [issue](https://github.com/scottcburgess/Bayesian-evolvability-NCII/issues). <br /> 
For questions about the paper, please contact [Scott Burgess](https://github.com/scottcburgess). 

