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

Additional packages needed to summarize and plot the posteriors of the model <br />
`modeest` <br />
`ggridges` <br />
`pander` <br />
`gridExtra` <br />

### Files <br />
`0_misc_funcs.R` miscellaneous functions used in summary and figure making files.  <br />
`1_load_and_format.R` load raw data and create formatted RDS files for use in models.  <br />

`2a_run_model_jags.R` script for running the main model in JAGS.  <br />

`Figure_#.R` scripts for making figures in manuscript. <br />
`posterior_summary.rmd` Rmarkdown file for summarizing posterior output from model runs.  <br />
`Table_1.rmd` Rmarkdown for making Table 1 in the manuscript.  <br />
`Supplementary Material 1.rmd` Rmarkdown file for creating Supplementary Material 1.  <br />

**Data** <br />
`Moccidentalis All Traits.csv` raw data file. <br />
`trunk_tail_data.rds` formatted data for Models I, II, and IV in RDS format.  <br />
`hatch_settle.rds` formatted data for Models III and IV in RDS format.  <br />

**Model_outputs** <br />
Contains the `*.rdata` file produced by `2a_run_model_jags.R`. The file size is too large to be stored on GitHub, so email us for a copy. <br />

**Figure_and_Tables** <br />
`Figure_#.pdf` figures in the paper  <br />
`*.html` Summary of the posteriors of the JAGS model  <br />
`Supplementary-Material-1.pdf` Output file  <br />
`Table_1.docx` Output file   <br />

### Contact
For questions about the Model code, please contact [Eric Archer](https://github.com/EricArcher). <br />
For questions about the Figure code, please contact [Scott Burgess](https://github.com/scottcburgess). <br />
For questions about the paper, please contact [Scott Burgess](https://github.com/scottcburgess).  <br />
For bug reports, please leave an [issue](https://github.com/scottcburgess/Bayesian-evolvability-NCII/issues). <br /> 

