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

### Files
**1_Data** <br />
`Moccidentalis All Traits.csv` raw data file. <br />
`Moccidentalis All Traits_updated.csv` has the dish ID added, where it was recorded.    
`head_tail_data.rds` formatted data for Model I and II in RDS format (see `?readRDS`).  <br />
`hatch_settle.rds` formatted data for Models III and IV in RDS format (see `?readRDS`).  <br />

**2_Model_R_scripts** <br />
`Model_##.R` scripts for running Models I - IV.  <br />
`0_misc_funcs.R` miscellaneous functions used in summary and figure making files.  <br />
`1_load_and_format.R` load raw data and create formatted RDS files for use in models.  <br />
`1_load_and_format_updated.R` (UNFINISHED) updated script to include the dish ID in `Moccidentalis All Traits_updated.csv`.   
`Model_I_MCMCglmm.R` is `Model_I.R` written in MCMCglmm for comparison.  
`Model_I_MCMCglmm_blocks_random.R` is `Model_I.R` written in MCMCglmm but with blocks as a random effect, for comparison.  
`Model_IV_MCMCglmm.R` is `Model_IV.R` re-imagined, based on reviewer comments, 
to fit trunk length, tail length, and probability of settlement as a single G-matrix.  

**3_Model_outputs** <br />
`Model_##.Rmd` Rmarkdown file for summarizing posterior output from model runs.  <br />
`Model_#_posterior_########_####.rdata` Image of all objects and posteriors produced by runs of Models I - IV.  <br />

**4_Figure_R_scripts** <br />
`Figure_#.R` scripts for making figures in manuscript.  <br />
`Results_Summary_and_Tables.R` scripts for making the summaries and tables in manuscript

**5_Figure_outputs** <br />
`Figure_#.pdf` figures in the paper


### Contact
For questions about the Model code, please contact [Eric Archer](https://github.com/EricArcher). <br />
For questions about the Figure code, please contact [Scott Burgess](https://github.com/scottcburgess). <br />
For questions about the paper, please contact [Scott Burgess](https://github.com/scottcburgess).  <br />
For bug reports, please leave an [issue](https://github.com/scottcburgess/Bayesian-evolvability-NCII/issues). <br /> 

