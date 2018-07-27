pkgname <- "metalr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('metalr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("forest_metalr")
### * forest_metalr

flush(stderr()); flush(stdout())

### Name: forest_metalr
### Title: Forest plot for likeliihood ratio based meta-analysis.
### Aliases: forest_metalr

### ** Examples

## Not run: 
##D data("statindata")  #statin potency and acute kidney injury dataset
##D # the metalr object
##D metalr_obj<-metalr_or(idata=statindata[,2:5],refval=0,num_iter=3000,increm=0.001,method = "random")
##D #forest plot of the metalr object
##D forest_metalr(metalr_obj)
## End(Not run)




cleanEx()
nameEx("ici.or")
### * ici.or

flush(stderr()); flush(stdout())

### Name: ici.or
### Title: 95% Intrinsic Confidence Interval (ICI) for Odds Ratio (OR) in
###   observational studies.
### Aliases: ici.or
### Keywords: ICIs Likelihood-ratio Meta-analysis

### ** Examples

## Not run: 
##D data("statindata") # statin potency and acute kidney injury data
##D ici.or(idata = statindata[1,2:5]) # ICI for study
## End(Not run)




cleanEx()
nameEx("ici.rr")
### * ici.rr

flush(stderr()); flush(stdout())

### Name: ici.rr
### Title: 95% Intrinsic confidence intervals for Rate Ratios (RR) in
###   epidemiological studies.
### Aliases: ici.rr
### Keywords: Likelihood-ratio Meta-analysis

### ** Examples

## Not run: 
##D # Clopidogrel vs Aspirin trial dataset
##D cases<-c(939,1021)
##D person_yrs<-c(17636,17519)
##D patients<-c(9599,9586)
##D ici.rr(cases, patients, person_yrs)
## End(Not run)




cleanEx()
nameEx("metalr_or")
### * metalr_or

flush(stderr()); flush(stdout())

### Name: metalr_or
### Title: Likelihood ratio meta-analysis for combining odds ratios in
###   fixed and random effects meta-analyses.
### Aliases: metalr_or
### Keywords: Likelihood-ratio Meta-analysis

### ** Examples

## Not run: 
##D # statin potency and acute kidney injury data
##D data("statindata")
##D metalr_or(idata=statindata[,2:5],refval=0,num_iter=3000,increm=0.001,method = "random")
## End(Not run)




cleanEx()
nameEx("metalr_rr")
### * metalr_rr

flush(stderr()); flush(stdout())

### Name: metalr_rr
### Title: Likelihood ratio meta-analysis for combining rate ratios in
###   fixed and random effects meta-analyses.
### Aliases: metalr_rr
### Keywords: Likelihood-ratio Meta-analysis

### ** Examples

## Not run: 
##D Random dataset
##D data("sample_metarr_data")
##D metalr_rr(idata=sample_metarr_data,refval=0,num_iter=3000,increm=0.001,method = "random")
## End(Not run)




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
