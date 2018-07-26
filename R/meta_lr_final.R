###  Functions for likelihood ratio meta-analysis R package "metalr"
###  to generate traditional and intrinsic confidence intervals for
###  odds ratio (OR) and rate ratio (RR)
##########################################################################


#' 95\% Intrinsic Confidence Interval (ICI) for Odds Ratio (OR) in observational studies.
#'
#' @description Calculates traditional and intrinsic confidence intervals for
#'  odds ratio from an observational study.
#'
#' @param idata Vector of length 4: cases for treatment A, controls for treatment A,
#'  cases for treatment B and control for treatment B.
#'
#' @return OR : MLE estimate of the odds ratio
#' @return llci : Lower 95\% traditional confidence limit
#' @return ulci : Upper 95\% traditional confidence limit
#' @return llici : Lower 95\% intrinsic confidence limit
#' @return ulici : Upper 95\% intrinsic confidence limit
#'
#' @references Dormuth, Colin R., Kristian B. Filion, and Robert W.
#' Platt. "Likelihood ratio meta-analysis: New motivation and approach
#'  for an old method." Contemporary clinical trials 47 (2016): 259-265.
#'
#' @examples
#' \dontrun{
#' data("statindata") # statin potency and acute kidney injury data
#' ici.or(idata = statindata[1,2:5]) # ICI for study
#' }
#'
#' @keywords Likelihood-ratio
#' @keywords ICIs
#' @keywords Meta-analysis
#'
#' @export
ici.or<-function(idata){
  #Data: case-controls pairs. eg. c(case1_trtA,ctrl_trtA,case2_trtB,ctrl_trtB)

  # Warning messages
  if(length(idata)<4){
    stop("The number of columns is less than 4. Please enter a vector of length 4. ")
  }
  if(length(idata)>4){
    warning("The number of columns is greater than 4. Only columns 1-4 of the imputed dataset are used. ")
  }

  tempdata<-as.numeric(idata)
  OR<-(tempdata[1]/tempdata[2])/(tempdata[3]/tempdata[4])
  log.OR<-log(OR)
  var.log.mle<-(1/tempdata[1])+(1/tempdata[2])+(1/tempdata[3])+(1/tempdata[4])

  ll<-exp(log.OR-1.96*(var.log.mle^(1/2)))
  ul<-exp(log.OR+1.96*(var.log.mle^(1/2)))
  ci<-round(c(ll,ul),3);#ci  #the confidence limits

  # the standard error of effect estimate
  log.ll<-log(ll)
  log.ul<-log(ul)
  SE<-(log.ul-log.ll)/(2*1.96)
  #var.log.mle=SE*SE; var.log.mle  # Just the same as the delta method

  # support for mle vs null
  supmle<-(log.OR^2)/(2*var.log.mle) #supmle

  # ratio of variance under the null and alternative
  vratio<-(2*var.log.mle)/(2*var.log.mle)

  # log of intrinsic confidence limits
  logHAL<-log.OR-sqrt(vratio*abs(log.OR^2-(supmle-2.99)*2*var.log.mle))
  logHAU<-(2*log.OR)-logHAL
  ici.ll<-exp(logHAL)
  ici.ul<-exp(logHAU)
  ORvec<-data.frame(OR=OR,ll.ci=ll,ul.ci=ul,ll.ici=ici.ll,ul.ici=ici.ul)
  return(ORvec)
}


#' 95\% Intrinsic confidence intervals for Rate Ratios (RR) in epidemiological studies.
#'
#' @description Calculates 95\% traditional confidence limits and 95\% intrinsic confidence
#' intervals for rate ratio from epidemiological studies.
#'
#' @param cases The number of individuals affected by the condition
#' @param patients The total number of individuals in the study
#' @param person_yrs The amount of time the patients were followed during the study
#'
#' @return RR MLE : estimate of the rates ratio
#' @return llci : Lower 95\% traditional confidence limit
#' @return ulci : Upper 95\% traditional confidence limit
#' @return llici : Lower 95\% intrinsic confidence limit
#' @return ulici : Upper 95\% intrinsic confidence limit
#'
#' @references Dormuth, Colin R., Kristian B. Filion, and Robert W.
#' Platt. "Likelihood ratio meta-analysis: New motivation and approach
#'  for an old method." Contemporary clinical trials 47 (2016): 259-265.
#'
#' @examples
#' \dontrun{
#' # Clopidogrel vs Aspirin trial dataset
#' cases<-c(939,1021)
#' person_yrs<-c(17636,17519)
#' patients<-c(9599,9586)
#' ici.rr(cases, patients, person_yrs)
#' }
#'
#' @keywords Likelihood-ratio
#' @keywords Meta-analysis
#'
#' @export
ici.rr<-function(cases,patients,person_yrs){
  #Data:vectors of group A and R respectively

  # Warning messages for incorrect entries...
  if(class(cases)!="numeric" | class(person_yrs)!="numeric" ){
    stop("cases and person-years contain non-numeric entries. Please enter numeric values.")
  }

  if(length(cases)==2 & length(person_yrs)==2){
    rate.ratio<-round((cases/person_yrs)[1]/(cases/person_yrs)[2],3)
    log.RR<- round(log(rate.ratio),3)
    var.log.mle<- round((1/cases[1])+(1/cases[2]),3)

    ll<-exp(log.RR-1.96*(var.log.mle^(1/2)))
    ul<-exp(log.RR+1.96*(var.log.mle^(1/2)))
    ci<-round(c(ll,ul),3);  #the confidence limits ci
    #1/exp(1.96)*100 # unit changes for a 95% CI is ~14% as opposed to 5% change expected

    # the standard error of effect estimate
    log.ll<-log(ll)
    log.ul<-log(ul)
    SE<-(log.ul-log.ll)/(2*1.96); SE
    #var.log.mle=SE*SE; var.log.mle  # Just the same as the delta method

    # support for mle vs null
    supmle<-(log.RR^2)/(2*var.log.mle);#supmle

    # ratio of variance under the null and alternative
    vratio<-(2*var.log.mle)/(2*var.log.mle)

    # log of intrinsic confidence limits
    logHAL<-log.RR-sqrt(vratio*abs(log.RR^2-(supmle-2.99)*2*var.log.mle))
    logHAU<-(2*log.RR)-logHAL
    ici.ll<-exp(logHAL)
    ici.ul<-exp(logHAU)
    RRvec<-data.frame(RR=rate.ratio,ll.ci=ll,ul.ci=ul,ll.ici=ici.ll,ul.ici=ici.ul)
    return(RRvec)
  }
  else{
    print("Function requires cases for treatment A&B and person years for A&B")
  }
}



### Defining the meta likelihood_ratio for odds ratio function
# Take a dataset of 4 columns. case and control for treatment A followed
# case and control for treatment B

#' Likelihood ratio meta-analysis for combining odds ratios in fixed
#' and random effects meta-analyses.
# Estimates 95\% traditional CIs and 95\% Intrinsic CIs of combined odds ratio.
#'
#' @param idata A dataframe of 4 columns for cases control pairs for treatments
#' @param refval The reference value for the log of the alternate hypothesis
#' @param num_iter The number of iterations or steps from the alternate hypothesis
#' @param increm The quantity of increments of the refval upto the number of iterations
#' @param method The meta-analytic method i.e. fixed or random effect method.
#'
#' @description Based on the method proposed by Dormuth et al, 2016,
#' the function estimates traditional 95\% confidence intervals and intrinsic
#' confidence intervals for combined effect estimates (OR) in meta-analysis.
#'  It also returns an estimate of heterogeneity accross studies as well as
#'  Isq statistics in random meta-analysis.
#'
#' @return Total_RE : A dataframe of total effect estimate from meta analysis,
#'  the 95\% CIs and intrinsic CIs.
#' @return  Tausq : Measure of heterogeneity between the studies used in the
#' meta-analysis in random effect meta-analysis.
#' @return Isq : The I^2 statistics
#' @return meta_results : Dataframe effect estimates from all the studies,
#'  the 95\% confidence limits and the 95\% intrinsic confidence limits.
#'
#'
#' @references Dormuth, Colin R., Kristian B. Filion, and Robert W.
#' Platt. "Likelihood ratio meta-analysis: New motivation and approach
#' for an old method." Contemporary clinical trials 47 (2016): 259-265.
#'
#' @examples
#' \dontrun{
#' # statin potency and acute kidney injury data
#' data("statindata")
#' metalr_or(idata=statindata[,2:5],refval=0,num_iter=3000,increm=0.001,method = "random")
#' }
#'
#' @keywords Likelihood-ratio
#' @keywords Meta-analysis
#'
#' @export
metalr_or<-function(idata,refval,num_iter,increm,method="random"){

  if(class(idata) != "data.frame"){
    stop("You did not supply the data as a dataframe. Please supply your data as a dataframe.")
  }
  if(ncol(idata) > 4){
    warning("The number of columns is greater than 4. Only columns 1-4 of the imputed dataset are used. ")
  }
  for (col in 1:4) {
    if(class(idata[,col])!="numeric")
      stop("Columns 1-4 contain non numeric entries. Please enter numeric values.")
  }

  if(class(refval) != "numeric" | length(refval)!= 1){
    stop("Please enter a single numeric entry for the 'refval' ")
  }

  if(class(num_iter) != "numeric" | length(num_iter)!= 1){
    stop("Please enter a single numeric entry for the number of iterations 'num_iter' ")
  }

  if(class(increm) != "numeric" | length(increm)!= 1){
    stop("Please enter a single numeric entry for the number of increments 'increm' ")
  }

  if (method%in%c("random","fixed")){
    num_study<-dim(idata)[1]
    # Weights  for meta-analysis
    tempdata<-idata
    tempdata$OR<-(tempdata[,1]/tempdata[,2])/(tempdata[,3]/tempdata[,4])
    tempdata$log.OR<-log(tempdata$OR)
    tempdata$var.log.mle<-(1/tempdata[,1])+(1/tempdata[,2])+(1/tempdata[,3])+(1/tempdata[,4])

    tempdata$wi<-1/tempdata$var.log.mle
    tempdata$wisq<-tempdata$wi^2
    tempdata$tisq<-tempdata$log.OR^2
    tempdata$witisq<-tempdata$wi*tempdata$tisq
    tempdata$witi<-tempdata$wi*tempdata$log.OR

    Q<-sum(tempdata$witisq)-((sum(tempdata$witi)^2)/sum(tempdata$wi))
    C<-sum(tempdata$wi)-(sum(tempdata$wisq)/sum(tempdata$wi))
    recs<-num_study

    # some modifications##
    tausq<-(Q-(recs-1))/C
    isq<-(Q-(recs-1))/Q

    # For the fixed effects components
    logLR.l<-rep(NA,num_iter)
    logLR.r<-rep(NA,num_iter)

    # For the random components
    rlogLR.l<-matrix(rep(NA,num_iter*(num_study+1)),nrow = num_iter)
    rlogLR.r<-matrix(rep(NA,num_iter*(num_study+1)),nrow = num_iter)


    LH1<-rep(NA,num_iter)
    RH1<-rep(NA,num_iter)

    log.OR<-tempdata$log.OR
    var.log.mle<-tempdata$var.log.mle

    for (i in 1:num_iter) {#search the neg and positive side of the ref H1
      loglr.l<-rep(NA,num_study)
      loglr.r<-rep(NA,num_study)
      rloglr.l<-rep(NA,num_study)
      rloglr.r<-rep(NA,num_study)
      for (s in 1:num_study) {
        lH1<-refval-(increm*i)
        loglr.l[s]<-((log.OR[s]^2)/(2*var.log.mle[s]))-((log.OR[s]-lH1)^2/(2*var.log.mle[s]))
        rloglr.l[s]<-(var.log.mle[s]/(var.log.mle[s]+tausq))*loglr.l[s]

        rH1<-refval+(increm*i)
        loglr.r[s]<-((log.OR[s]^2)/(2*var.log.mle[s]))-((log.OR[s]-rH1)^2/(2*var.log.mle[s]))
        rloglr.r[s]<-(var.log.mle[s]/(var.log.mle[s]+tausq))*loglr.r[s]
      }
      logLR.l[i]<-sum(loglr.l)
      logLR.r[i]<-sum(loglr.r)
      rlogLR.l[i,]<-c(loglr.l,sum(rloglr.l))
      rlogLR.r[i,]<-c(loglr.r,sum(rloglr.r)) # we will search over these possible values
      LH1[i]<-lH1
      RH1[i]<-rH1 # lets append the values of H_alt
    }
    colnames(rlogLR.l)<-c(1:num_study,"total_relk")
    colnames(rlogLR.r)<-c(1:num_study,"total_relk")

    # check the function here
    res.lk<-data.frame(H1=c(LH1,RH1), rbind(rlogLR.l,rlogLR.r),
                       "total_felk"=c(logLR.l,logLR.r))
    ###########
    # Fixed Effect LR_mle values
    felrmax<-max(res.lk$total_felk)
    feH1<-res.lk$H1[which.max(res.lk$total_felk)]

    # search for lower limit of 95% ICI and CI for the fixed effect
    llicidata<-res.lk[res.lk$total_felk<=(felrmax-2.99)&res.lk$H1<feH1,]
    ll_log_ici<-max(llicidata$total_felk)
    ll_ici_H1<-exp(llicidata$H1[which.max(llicidata$total_felk)])

    llcidata<-res.lk[res.lk$total_felk<=(felrmax-1.96)&res.lk$H1<feH1,]
    ll_log_ci<-max(llcidata$total_felk)
    ll_ci_H1<-exp(llcidata$H1[which.max(llcidata$total_felk)])


    # search for upper limit of 95% ICI and CI
    ulicidata<-res.lk[res.lk$total_felk<=(felrmax-2.99)&res.lk$H1>feH1,]
    ul_log_ici<-max(ulicidata$total_felk)
    ul_ici_H1<-exp(ulicidata$H1[which.max(ulicidata$total_felk)])

    ulcidata<-res.lk[res.lk$total_felk<=(felrmax-1.96)&res.lk$H1>feH1,]
    ul_log_ci<-max(ulcidata$total_felk)
    ul_ci_H1<-exp(ulcidata$H1[which.max(ulcidata$total_felk)])

    ## the ci and ici limits for FE
    fixed_ci<-c(exp(feH1),ll_ci_H1,ul_ci_H1)
    fixed_ici<-c(exp(feH1),ll_ici_H1,ul_ici_H1)


    #############
    # Random effect lr_mle values | relrmax
    relrmax<-max(res.lk$total_relk)
    reH1<-res.lk$H1[which.max(res.lk$total_relk)]

    # search for lower limit of 95% ICI and CI for the random effect
    rellicidata<-res.lk[res.lk$total_relk<=(relrmax-2.99)&res.lk$H1<reH1,]
    rell_log_ici<-max(rellicidata$total_relk)
    rell_ici_H1<-exp(rellicidata$H1[which.max(rellicidata$total_relk)])

    rellcidata<-res.lk[res.lk$total_relk<=(relrmax-1.96)&res.lk$H1<reH1,]
    rell_log_ci<-max(rellcidata$total_relk)
    rell_ci_H1<-exp(rellcidata$H1[which.max(rellcidata$total_relk)])

    # search for upper limit of 95% ICI and CI
    reulicidata<-res.lk[res.lk$total_relk<=(relrmax-2.99)&res.lk$H1>reH1,]
    reul_log_ici<-max(reulicidata$total_relk)
    reul_ici_H1<-exp(reulicidata$H1[which.max(reulicidata$total_relk)])

    reulcidata<-res.lk[res.lk$total_relk<=(relrmax-1.96)&res.lk$H1>reH1,]
    reul_log_ci<-max(reulcidata$total_relk)
    reul_ci_H1<-exp(reulcidata$H1[which.max(reulcidata$total_relk)])

    ## the ci and ici limits for Random effects
    rand_ci<-c(exp(reH1),rell_ci_H1,reul_ci_H1)
    rand_ici<-c(exp(reH1),rell_ici_H1,reul_ici_H1)


    # lets compile all the individual estimates of ci and intrinsic cis

    out_effects<-data.frame(study=c(1:num_study,"Total"),MLE=rep(NA,num_study+1),
                            llci=rep(NA,num_study+1),ulci=rep(NA,num_study+1),
                            llici=rep(NA,num_study+1),ulici=rep(NA,num_study+1) )
    for (k in 1:num_study) {
      out_effects[k,2:6]<-ici.or(idata = idata[k,] )
    }

    if (method=="fixed") {
      out_effects[num_study+1,2:6]<-c(fixed_ci,fixed_ici[2:3])
      return(list(Total_FE=data.frame(limit=c("MLE","L95%","U95%"),CI=fixed_ci,ICI=fixed_ici),
                  meta_result=out_effects))
    }
    if (method=="random") {
      out_effects[num_study+1,2:6]<-c(rand_ci,rand_ici[2:3])
      return(list(Total_RE=data.frame(limit=c("MLE","L95%","U95%"),CI=rand_ci,ICI=rand_ici),
                  TauSq=tausq, Isq=isq, meta_result=out_effects))
    }
  }
  else{
    print("Please select a 'fixed' or 'random' effect method")
  }
}




##########################################################
# Define the meta likelihood ratio for Rate-Ratio
# Takes a dataset/dataframe with each row representing a study
# six columns: case.trt.A, case.trt.B, pers.time.trt.A, pers.time.trt.B,
# patients.A, patients.B... see sampledataset

#' Likelihood ratio meta-analysis for combining rate ratios in fixed
#' and random effects meta-analyses.
# Estimates 95\% traditional CIs and 95\% Intrinsic CIs of combined odds ratio.
#'
#' @param idata A dataframe of atleast 4 columns of: cases for treatment A, cases for
#' treatment B, person time for treatment A and person time for treatment B.
#' @param refval The reference value for the log of the alternate hypothesis
#' @param num_iter The number of iterations or steps from the alternate hypothesis
#' @param increm The quantity of increments of the refval upto the number of iterations
#' @param method The meta-analytic method i.e. fixed or random effect method.
#'
#' @description Based on the method proposed by Dormuth et al, 2016,
#' the function estimates traditional 95\% confidence intervals and intrinsic
#' confidence intervals for combined effect estimates (RR) in meta-analysis.
#'  It also returns an estimate of heterogeneity accross studies as well as
#'  Isq statistics in random meta-analysis.
#'
#' @return Total_RE : A dataframe of total effect estimate from meta analysis,
#'  the 95\% CIs and intrinsic CIs.
#' @return  Tausq : Measure of heterogeneity between the studies used in the
#' meta-analysis in random effect meta-analysis.
#' @return Isq : The I^2 statistics
#' @return meta_results : Dataframe effect estimates from all the studies,
#'  the 95\% confidence limits and the 95\% intrinsic confidence limits.
#'
#'
#' @references Dormuth, Colin R., Kristian B. Filion, and Robert W.
#' Platt. "Likelihood ratio meta-analysis: New motivation and approach
#' for an old method." Contemporary clinical trials 47 (2016): 259-265.
#'
#' @examples
#' \dontrun{
#' Random dataset
#' data("sample_metarr_data")
#' metalr_rr(idata=sample_metarr_data,refval=0,num_iter=3000,increm=0.001,method = "random")
#' }
#'
#' @keywords Likelihood-ratio
#' @keywords Meta-analysis
#'
#' @export
metalr_rr<-function(idata,refval,num_iter,increm,method="random"){

  # error and warning messages
  if(class(idata) != "data.frame"){
    stop("You did not supply the data as a dataframe. Please supply your data as a dataframe.")
  }
  if(ncol(idata) > 4){
    warning("The number of columns is greater than 4. Only columns 1-4 of the imputed dataset are used. ")
  }
  for (col in 1:4) {
    if(class(idata[,col])!="numeric")
      stop("Columns 1-4 contain non numeric entries. Please enter numeric values.")
  }

  if(class(refval) != "numeric" | length(refval)!= 1){
    stop("Please enter a single numeric entry for the 'refval' ")
  }

  if(class(num_iter) != "numeric" | length(num_iter)!= 1){
    stop("Please enter a single numeric entry for the number of iterations 'num_iter' ")
  }

  if(class(increm) != "numeric" | length(increm)!= 1){
    stop("Please enter a single numeric entry for the number of increments 'increm' ")
  }


  if (method%in%c("random","fixed")){
    num_study<-dim(idata)[1]

    # Weights  for meta-analysis
    tempdata<-idata
    tempdata$RR<-(tempdata[,1]/tempdata[,3])/(tempdata[,2]/tempdata[,4])
    tempdata$log.RR<-log(tempdata$RR)
    tempdata$var.log.mle<-(1/tempdata[,1])+(1/tempdata[,2])

    tempdata$wi<-1/tempdata$var.log.mle
    tempdata$wisq<-tempdata$wi^2
    tempdata$tisq<-tempdata$log.RR^2
    tempdata$witisq<-tempdata$wi*tempdata$tisq
    tempdata$witi<-tempdata$wi*tempdata$log.RR

    Q<-sum(tempdata$witisq)-((sum(tempdata$witi)^2)/sum(tempdata$wi))
    C<-sum(tempdata$wi)-(sum(tempdata$wisq)/sum(tempdata$wi))
    recs<-num_study

    # some modifications##
    tausq<-(Q-(recs-1))/C
    isq<-(Q-(recs-1))/Q

    # For the fixed effects components
    logLR.l<-rep(NA,num_iter)
    logLR.r<-rep(NA,num_iter)

    # For the random components
    rlogLR.l<-matrix(rep(NA,num_iter*(num_study+1)),nrow = num_iter)
    rlogLR.r<-matrix(rep(NA,num_iter*(num_study+1)),nrow = num_iter)


    LH1<-rep(NA,num_iter)
    RH1<-rep(NA,num_iter)

    log.RR<-tempdata$log.RR
    var.log.mle<-tempdata$var.log.mle

    for (i in 1:num_iter) {#search the neg and positive side of the ref H1
      loglr.l<-rep(NA,num_study)
      loglr.r<-rep(NA,num_study)
      rloglr.l<-rep(NA,num_study)
      rloglr.r<-rep(NA,num_study)
      for (s in 1:num_study) {
        lH1<-refval-(increm*i)
        loglr.l[s]<-((log.RR[s]^2)/(2*var.log.mle[s]))-((log.RR[s]-lH1)^2/(2*var.log.mle[s]))
        rloglr.l[s]<-(var.log.mle[s]/(var.log.mle[s]+tausq))*loglr.l[s]

        rH1<-refval+(increm*i)
        loglr.r[s]<-((log.RR[s]^2)/(2*var.log.mle[s]))-((log.RR[s]-rH1)^2/(2*var.log.mle[s]))
        rloglr.r[s]<-(var.log.mle[s]/(var.log.mle[s]+tausq))*loglr.r[s]
      }
      logLR.l[i]<-sum(loglr.l)
      logLR.r[i]<-sum(loglr.r)
      rlogLR.l[i,]<-c(loglr.l,sum(rloglr.l))
      rlogLR.r[i,]<-c(loglr.r,sum(rloglr.r)) # we will search over these possible values
      LH1[i]<-lH1
      RH1[i]<-rH1 # lets append the values of H_alt
    }
    colnames(rlogLR.l)<-c(1:num_study,"total_relk")
    colnames(rlogLR.r)<-c(1:num_study,"total_relk")

    # check the function here
    res.lk<-data.frame(H1=c(LH1,RH1), rbind(rlogLR.l,rlogLR.r),
                       "total_felk"=c(logLR.l,logLR.r))
    ###########
    # Fixed Effect LR_mle values
    felrmax<-max(res.lk$total_felk)
    feH1<-res.lk$H1[which.max(res.lk$total_felk)]

    # search for lower limit of 95% ICI and CI for the fixed effect
    llicidata<-res.lk[res.lk$total_felk<=(felrmax-2.99)&res.lk$H1<feH1,]
    ll_log_ici<-max(llicidata$total_felk)
    ll_ici_H1<-exp(llicidata$H1[which.max(llicidata$total_felk)])

    llcidata<-res.lk[res.lk$total_felk<=(felrmax-1.96)&res.lk$H1<feH1,]
    ll_log_ci<-max(llcidata$total_felk)
    ll_ci_H1<-exp(llcidata$H1[which.max(llcidata$total_felk)])


    # search for upper limit of 95% ICI and CI
    ulicidata<-res.lk[res.lk$total_felk<=(felrmax-2.99)&res.lk$H1>feH1,]
    ul_log_ici<-max(ulicidata$total_felk)
    ul_ici_H1<-exp(ulicidata$H1[which.max(ulicidata$total_felk)])

    ulcidata<-res.lk[res.lk$total_felk<=(felrmax-1.96)&res.lk$H1>feH1,]
    ul_log_ci<-max(ulcidata$total_felk)
    ul_ci_H1<-exp(ulcidata$H1[which.max(ulcidata$total_felk)])

    ## the ci and ici limits for FE
    fixed_ci<-c(exp(feH1),ll_ci_H1,ul_ci_H1)
    fixed_ici<-c(exp(feH1),ll_ici_H1,ul_ici_H1)

    #############
    # Random effect lr_mle values | relrmax
    relrmax<-max(res.lk$total_relk)
    reH1<-res.lk$H1[which.max(res.lk$total_relk)]

    # search for lower limit of 95% ICI and CI for the random effect
    rellicidata<-res.lk[res.lk$total_relk<=(relrmax-2.99)&res.lk$H1<reH1,]
    rell_log_ici<-max(rellicidata$total_relk)
    rell_ici_H1<-exp(rellicidata$H1[which.max(rellicidata$total_relk)])

    rellcidata<-res.lk[res.lk$total_relk<=(relrmax-1.96)&res.lk$H1<reH1,]
    rell_log_ci<-max(rellcidata$total_relk)
    rell_ci_H1<-exp(rellcidata$H1[which.max(rellcidata$total_relk)])

    # search for upper limit of 95% ICI and CI
    reulicidata<-res.lk[res.lk$total_relk<=(relrmax-2.99)&res.lk$H1>reH1,]
    reul_log_ici<-max(reulicidata$total_relk)
    reul_ici_H1<-exp(reulicidata$H1[which.max(reulicidata$total_relk)])

    reulcidata<-res.lk[res.lk$total_relk<=(relrmax-1.96)&res.lk$H1>reH1,]
    reul_log_ci<-max(reulcidata$total_relk)
    reul_ci_H1<-exp(reulcidata$H1[which.max(reulcidata$total_relk)])

    ## the ci and ici limits for Random effects
    rand_ci<-c(exp(reH1),rell_ci_H1,reul_ci_H1)
    rand_ici<-c(exp(reH1),rell_ici_H1,reul_ici_H1)

    # lets compile all the individual estimates of ci and intrinsic cis
    cases = cbind(tempdata[,1],tempdata[,2])
    person_yrs = cbind(tempdata[,3],tempdata[,4])

    out_effects<-data.frame(study=c(1:num_study,"Total"),MLE=rep(NA,num_study+1),
                            llci=rep(NA,num_study+1),ulci=rep(NA,num_study+1),
                            llici=rep(NA,num_study+1),ulici=rep(NA,num_study+1) )
    for (k in 1:num_study) {
      out_effects[k,2:6]<-ici.rr(cases = cases[k,],person_yrs = person_yrs[k,])
    }


    if (method=="fixed") {
      out_effects[num_study+1,2:6]<-c(fixed_ci,ll_ici_H1,ul_ici_H1)
      return(list(Total_FE=data.frame(limit=c("MLE","L95%","U95%"),CI=fixed_ci,ICI=fixed_ici),
                  meta_result=out_effects))
    }
    if (method=="random") {
      out_effects[num_study+1,2:6]<-c(rand_ci,rand_ici[2:3])
      return(list(Total_RE=data.frame(limit=c("MLE","L95%","U95%"),CI=rand_ci,ICI=rand_ici),
                  TauSq=tausq, Isq=isq, meta_result=out_effects))
    }
  }
  else{
    print("select a 'fixed' or 'random' effect method!")
  }
}


##############################################
# Adding a forrest plot option to the metalr object


#' Forest plot for likeliihood ratio based meta-analysis.
# showing 95\% traditional CIs and 95\% Intrinsic CIs.
#'
#' @description The function plots confidence limits of traditional 95\% CIs
#' and 95\% ICIs for the studies included in the meta-analysis as well as confidence bars
#' associated with the overall effect estimates.
#'
#' @param metalr_obj An abject from the metalr functions metalr_or()
#'  or metalr_rr(). The metalr object is a list of results computed
#'   by the metalr functions which includes a dataframe of mle of the
#'   effect estimates and their corresponding 95\% CIs and ICIs.
#'   See the example below.
#'
#' @return Returns a forest plot of the 95\% CIs and 95\% ICIs.

#' @examples
#' \dontrun{
#' data("statindata")  #statin potency and acute kidney injury dataset
#' # the metalr object
#' metalr_obj<-metalr_or(idata=statindata[,2:5],refval=0,num_iter=3000,increm=0.001,method = "random")
#' #forest plot of the metalr object
#'forest_metalr(metalr_obj)
#' }
#'
#' @export
forest_metalr<-function(metalr_obj){
  metalr_obj_conf<-as.data.frame(metalr_obj$meta_result,row.names = NULL)
  row_names<-as.character(metalr_obj_conf$study)
  mle<-metalr_obj_conf$MLE
  lowerci<-metalr_obj_conf$llci
  upperci<-metalr_obj_conf$ulci
  lowerici<-metalr_obj_conf$llici
  upperici<-metalr_obj_conf$ulici
  metalr_plot<-forestplot::forestplot(labeltext=row_names,
                            mean=cbind(mle,mle),
                            lower=cbind(lowerci,lowerici),
                            upper=cbind(upperci,upperici),
                            txt_gp=forestplot::fpTxtGp(label=grid::gpar(cex=1.0),
                                           ticks=grid::gpar(cex=c(.9)),
                                           xlab= grid::gpar(cex = 0.8),
                                           title=grid::gpar(cex = 0.8)),
                            col=forestplot::fpColors(box=c("darkred","blue"), lines=c("red","blue"), zero = "gray50"),
                            zero=1, cex=0.9, lineheight = "auto", boxsize=0.2, colgap=grid::unit(3,"mm"),
                            lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2,lty.ci = c(1, 2),
                            xlab = "Effect Estimate",  clip =c(-1.125, 3.75),
                            title = "LRMA meta-analysis. \n Solid bars denote 95% CIs and dashed bars denote 95% ICIs")
  return(metalr_plot)
}

