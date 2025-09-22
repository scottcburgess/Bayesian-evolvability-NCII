# Bayesian-evolvability-NCII

Bayesian analysis of *Molgula occidentalis* (ascidian) larval traits from a North Carolina II (NCII) breeding design. <br />

### Reference
Powell JA, Archer FI, Burgess SC. In Prep. Evolvability of marine dispersal traits related to larval swimming and reproductive strategy.

### R Package Dependencies
`rjags` and local installation of [JAGS](https://mcmc-jags.sourceforge.io/) <br />
`evolvability` <br />
`QGglmm` <br />
`swfscMisc` <br />
`abind` <br />
`tidyverse` <br />

### Files <br />
`0_misc_funcs.R` miscellaneous functions used in summary and figure making files.  <br />
`1_load_and_format.R` load raw data and create formatted RDS files for use in models.  <br />

`Model_##.R` scripts for running Models I - IV.  <br />

`Figure_#.R` scripts for making figures in manuscript. <br />
`Table_1.rmd` Rmarkdown for making Table 1 in the manuscript.  <br />

**Data** <br />
`Moccidentalis All Traits.csv` raw data file. <br />
`trunk_tail_data.rds` formatted data for Models I, II, and IV in RDS format.  <br />
`hatch_settle.rds` formatted data for Models III and IV in RDS format.  <br />

**Model_outputs** <br />
`Model_##_summary.Rmd` Rmarkdown file for summarizing posterior output from model runs.  <br />

**Figure_outputs** <br />
`Figure_#.pdf` figures in the paper  <br />


### Contact
For questions about the Model code, please contact [Eric Archer](https://github.com/EricArcher). <br />
For questions about the Figure code, please contact [Scott Burgess](https://github.com/scottcburgess). <br />
For questions about the paper, please contact [Scott Burgess](https://github.com/scottcburgess).  <br />
For bug reports, please leave an [issue](https://github.com/scottcburgess/Bayesian-evolvability-NCII/issues). <br /> 

