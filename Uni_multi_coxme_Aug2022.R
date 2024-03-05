Compare_models <- function(baseModel, coxphModel,coxmeModel) {
  if (!class(coxphModel) == "coxph" | !class(coxmeModel) == "coxme") {
    stop("Wrong models")
  }
  
  anova_compare <- anova(baseModel, coxphModel, coxmeModel)
  print(anova_compare)
  
  ## Degrees of freedom
  coxphDf <- sum(!is.na(coef(coxphModel)))
  coxmeDf <- coxmeModel$df
  names(coxmeDf) <- c("Integrated","Penalized")
  
  ## DF differnces
  dfDiff <- coxmeDf - coxphDf
  ## Log likelihoods
  coxphLogLik <- coxphModel$loglik[2]  ## Log likelihood
  coxmeLogLik <- coxmeModel$loglik + c(0, 0, coxmeModel$penalty)
  ## Chi-squared value comparing the penalized model to the null model
  
  Loglikcompareph <-  -2 * (baseModel$loglik - logLik(coxphModel))
  Loglikcompareme <-  -2 * (coxmeLogLik["NULL"] - coxmeLogLik["Penalized"])
  
  outdf0 <- data.frame(c(Loglikcompareph,Loglikcompareme))
  rownames(outdf0) <- c("coxph to NULL model", "coxme to NULL model")
  colnames(outdf0) <- "Chi square value difference" 
  print(outdf0)
  
  ##______TESTING FOR RANDOM EFFECTS_____##
  ## -2 logLik difference
  logLikNeg2Diff <- c(-2 * (coxphLogLik - coxmeLogLik["Integrated"]),
                      -2 * (coxphLogLik - coxmeLogLik["Penalized"]))
  
  ## p-values
  pVals <- pchisq(q = logLikNeg2Diff, df = dfDiff, lower.tail = FALSE)
  
  ## Combine
  outDf <- data.frame(dfDiff, logLikNeg2Diff, pVals)
  colnames(outDf) <- c("Degrees of freedom difference", "-2logLik diff", "p-values")
  print(outDf)
  #formattable::formattable(anova_compare, align =c("l","c","c","c", "c","c","c", "c"))
}

filterCV <- function(data, PATID, names){
  cvlist <- data.frame()
  k = 1
  for (i in unique(data$PATID)) {
    my_data_subset <- data %>% filter(data$PATID==i)
    cv <- sd(my_data_subset[[names]]) / mean(my_data_subset[[names]]) * 100
    cvlist[k,1] <- i
    cvlist[k,2] <- cv
    k = k+1
  }
  cvlist <- as.data.frame(cvlist) %>% `colnames<-`(c("PATID", "CV"))
  cvlist$CV <- cvlist$CV %>%  replace_na(0)
  return(cvlist)
  
}

RunCoxPH <- function(data, names, Time = "OS", Event="OS.01"){
  f1 = as.formula(paste("Surv(", Time, ",", Event,") ~", paste(names, collapse="+")))
  skip_to_next <- FALSE
  tryCatch(coxph(f1, data=data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  md <- coxph(f1, data = data)
  obj <- univ_cox_ph(md)
  return(obj)
}

RunCoxME <- function(data, names, Time = "OS", Event="OS.01"){
  f1 = as.formula(paste("Surv(", Time, ",", Event,") ~", paste(names, "+ (1|PATID)")))
  skip_to_next <- FALSE
  tryCatch(coxme(f1, data=data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  md <- coxme(f1, data)
  obj <- univ_cox_me(md)
  
  return(obj)
}

## Univariate coxME model summarization 
univ_cox_me <- function(xx){
  beta <- xx$coefficients #$fixed is not needed
  hr <- exp(beta)
  z1 <- qnorm((1 + 0.95)/2, 0, 1)
  nvar <- length(beta)
  nfrail <- nrow(xx$var) - nvar
  se <- sqrt(diag(xx$var)[nfrail + 1:nvar])
  z <- round(beta/se, 2) 
  hr_low <- exp(beta - z1 * se)
  hr_high <- exp(beta + z1 * se)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  cc <- cox.zph(xx)
  tmp <- cc$table
  check_pval <- tmp[1,3]
  obj <- data.frame(cbind(nvar,nfrail,beta,hr,hr_low,hr_high,se,z,p,check_pval))
  return(obj)
}

## Univariate cox ph model summarization 
univ_cox_ph <- function(xx){
  aa <- summary(xx)
  tmp1 <- aa$conf.int
  tmp2 <- aa$coefficients
  ntot <- aa$n
  hr = tmp1[1,1]
  hr_low <- tmp1[1,3]
  hr_high <- tmp1[1,4]
  coeff <- tmp2[1,1]
  pval <- tmp2[1,5] ## for hazard ratio
  cc <- cox.zph(xx)
  tmp <- cc$table
  varnames <- rownames(tmp);
  p_check <- tmp[1,3]
  vname <- varnames[1]
  obj <- data.frame(cbind(vname,ntot,coeff, hr, hr_low, hr_high, pval, p_check))
  return(obj)
}

## Multivariate cox ph model summarization 
multi_cox_ph <- function(xx){
  aa <- summary(xx)
  tmp1 <- aa$conf.int
  tmp2 <- aa$coefficients
  ntot <- aa$n
  hr = tmp1[,1]
  hr_low <- tmp1[,3]
  hr_high <- tmp1[,4]
  coeff <- tmp2[,1]
  pval <- tmp2[,5] ## for hazard ratio
  cc <- cox.zph(xx)
  tmp <- cc$table
  varnames <- rownames(tmp);
  p_check <- tmp[,3]
  vname <- varnames[1]
  obj <- data.frame(cbind(vname,ntot,coeff, hr, hr_low, hr_high, pval, p_check))
  return(obj)
}

# name == name to be given on files 
# s = starting value of the omic dataset 
# e = ending value of the omic dataset 
# formula = formula in the main survival == Surv(time, event)~
# outpufilepath = Where the images are to be stored 
# PATID = Patient ID variable name
# SegmentLabel = Segment type / AOI type/ Cell type column header
# Time = Survival time 
# Event = Survival event at the end of survival time 
# printforest = Binary if forestplot is to be printed 
# multivariate = selection of which univariate model is to be chosen to run multivariate further 

Uni_multi_coxme <- function(data, name, s, e, uniquefactortype =NULL,significance, outputfilepath, PATID = "PATID", SegmentLabel="SegmentLabel", Time = "OS", Event="OS.01", printforest = FALSE, multivariatetype = c("CoxME", "CoxPH", "Aggregated Mean", "Aggregated Median"), filterCV = FALSE, CVcutoff = 30){
  
  vl <- names(data[,s:e])
  len <- length(vl)
  final_uni_coxph_table <- data.frame()
  final_uni_coxme_table <- data.frame()
  final_uni_aggregated_mean_table <- data.frame()
  final_uni_aggregated_median_table <- data.frame()
  
  
  # This chunk calculates coxph
  for (j in 1:len){
    if (filterCV == TRUE) {
      #print(vl[j])
      cvlist <- filterCV(data=data, PATID = "PATID", names=vl[j])
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data2 <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
      
      output <- RunCoxPH(data = data2, names=vl[j], Time = Time, Event=Event)
      final_uni_coxph_table <- final_uni_coxph_table %>% rbind(output)
      
    } else{
      output <- RunCoxPH(data = data, names=vl[j], Time = Time, Event= Event)
      final_uni_coxph_table <- final_uni_coxph_table %>% rbind(output)
    }
  }
  
  colnames(final_uni_coxph_table) <- c("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
  write.csv(final_uni_coxph_table, paste0("Cox regression outputs/",name,"_coxPH.csv"), append = TRUE)
  
  # This chunk calculates coxme
  for (j in 1:len){
    names = vl[j]
    if (filterCV == TRUE) {
      
      cvlist <- filterCV(data=data,  PATID = "PATID", names=vl[j])
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data2 <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
      f1 = as.formula(paste("Surv(", Time, ",", Event,") ~", paste(vl[j], "+ (1|PATID)")))
      skip_to_next <- FALSE
      tryCatch(coxme(f1, data=data2), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { next }
      md <- coxme(f1, data2)
      obj <- univ_cox_me(md)
      final_uni_coxme_table <- final_uni_coxme_table %>% rbind(obj)
      
    } else {
      
      f1 = as.formula(paste("Surv(", Time, ",", Event,") ~", paste(vl[j], "+ (1|PATID)")))
      skip_to_next <- FALSE
      tryCatch(coxme(f1, data=data), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { next }
      md <- coxme(f1, data)
      obj <- univ_cox_me(md)
      final_uni_coxme_table <- final_uni_coxme_table %>% rbind(obj)
    }
  }
  colnames(final_uni_coxme_table) <- c("nvar", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "SE", "Z", "Pvalue","Pvalue_check")
  final_uni_coxme_table <- final_uni_coxme_table %>%rownames_to_column("Biomarkers") 
  write.csv(final_uni_coxme_table, paste0("Cox regression outputs/",name,"_coxME.csv"), append = TRUE)
  
  final_uni_coxme_table <- final_uni_coxme_table %>% select("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
  
  ## This chunk calculates Aggregated Mean coxph
  
  for (j in 1:len){
    names = vl[j]
    
    if (filterCV == TRUE) {
      
      cvlist <- filterCV(data=data,  PATID = "PATID", names=names)
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data2 <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
      
      if (is.null(uniquefactortype)) {
        shortdf <-  data2 %>% dplyr::select(PATID,Time, Event, vl[j]) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, mean)
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_mean_table <- final_uni_aggregated_mean_table %>% rbind(output)
        # print(vl[j])
        
      } else {
        shortdf <-  data2 %>% dplyr::select(PATID, paste0(uniquefactortype),Time, Event,  vl[j]) %>% dplyr::group_by(PATID, paste0(uniquefactortype)) %>% summarise_if(is.numeric, mean)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_mean_table <- final_uni_aggregated_mean_table %>% rbind(output)
        
      }
    } else {
      if (is.null(uniquefactortype)) {
        shortdf <-  data %>% dplyr::select(PATID, Time, Event, vl[j]) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, mean)
        
        output <- RunCoxPH(data = shortdf, names=vl[j], Time = "OS", Event="OS.01")
        final_uni_aggregated_mean_table <- final_uni_aggregated_mean_table %>% rbind(output)
        
      } else {
        shortdf <-  data %>% dplyr::select(PATID,paste0(uniquefactortype), Time, Event, vl[j]) %>% dplyr::group_by(PATID, paste0(uniquefactortype))  %>% summarise_if(is.numeric, mean)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_mean_table <- final_uni_aggregated_mean_table %>% rbind(output)
        
      }
    }
  }
  
  colnames(final_uni_aggregated_mean_table) <- c("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
  write.csv(final_uni_aggregated_mean_table, paste0("Cox regression outputs/",name,"_A-Mean.csv"), append = TRUE)
  
  ## This chunk calculates Aggregated Median coxph
  for (j in 1:len){
    names = vl[j]
    
    if (filterCV == TRUE) {
      
      cvlist <- filterCV(data=data, PATID = "PATID", names=names)
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data2 <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
      
      if (is.null(uniquefactortype)) {
        shortdf <-  data2 %>% dplyr::select(PATID,Time, Event, vl[j]) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, median)
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_median_table <- final_uni_aggregated_median_table %>% rbind(output)
        # print(vl[j]) 
        
      } else {
        shortdf <-  data2 %>% dplyr::select(PATID,paste0(uniquefactortype),Time, Event,  vl[j]) %>% dplyr::group_by(PATID, paste0(uniquefactortype)) %>% summarise_if(is.numeric, median)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_median_table <- final_uni_aggregated_median_table %>% rbind(output)
        
      }
    } else {
      if (is.null(uniquefactortype)) {
        shortdf <-  data %>% dplyr::select(PATID, Time, Event, vl[j]) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, median)
        
        output <- RunCoxPH(data = shortdf, names=vl[j], Time = "OS", Event="OS.01")
        final_uni_aggregated_median_table <- final_uni_aggregated_median_table %>% rbind(output)
        
      } else {
        shortdf <-  data %>% dplyr::select(PATID,paste0(uniquefactortype), Time, Event, vl[j]) %>% dplyr::group_by(PATID, paste0(uniquefactortype))  %>% summarise_if(is.numeric, median)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_median_table <- final_uni_aggregated_median_table %>% rbind(output)
        
      }
    }
  }
  
  colnames(final_uni_aggregated_median_table) <- c("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
  write.csv(final_uni_aggregated_median_table, paste0("Cox regression outputs/",name,"_A-Median.csv"), append = TRUE)
  
  ## Plotting different tables
  
  total_table <- inner_join(final_uni_aggregated_median_table %>% `colnames<-`(c("Biomarkers", "n (A median)", "Beta (A median)","HR (A median)", "95% Lower HR (A median)", "95% Upper HR (A median)", "Pvalue (A median)","Pvalue_check (A median)")), final_uni_aggregated_mean_table %>% `colnames<-`(c("Biomarkers", "n (A mean)", "Beta (A mean)","HR (A mean)", "95% Lower HR (A mean)", "95% Upper HR (A mean)", "Pvalue (A mean)","Pvalue_check (A mean)")),  by = "Biomarkers")
  
  total_table <- inner_join(total_table, final_uni_coxph_table %>% `colnames<-`(c("Biomarkers", "n (PH)", "Beta (PH)","HR (PH)", "95% Lower HR (PH)", "95% Upper HR (PH)", "Pvalue (PH)","Pvalue_check (PH)")), by = "Biomarkers")
  
  total_table <- inner_join( total_table,  final_uni_coxme_table %>% `colnames<-`(c("Biomarkers", "n (CM)", "Beta (CM)","HR (CM)", "95% Lower HR (CM)", "95% Upper HR (CM)", "Pvalue (CM)","Pvalue_check (CM)")),  by = "Biomarkers")
  
  ## Creating shortlists
  shortlist0 <- final_uni_aggregated_median_table %>% select(n, Biomarkers, HR, Pvalue, Pvalue_check)
  shortlist1 <- final_uni_coxph_table %>% select(n, Biomarkers, HR, Pvalue, Pvalue_check)
  shortlist2 <- final_uni_coxme_table %>% select(n, Biomarkers, HR, Pvalue, Pvalue_check)
  shortlist3 <- final_uni_aggregated_mean_table %>% select(n, Biomarkers, HR, Pvalue, Pvalue_check)
  
  compare_table0 <- full_join(shortlist0,shortlist3, by="Biomarkers")
  colnames(compare_table0) <-  c( "n_med", "Biomarkers", "HR_med", "P_med", "P_chk_med", "n_mean", "HR_mean", "P_mean", "P_chk_mean")
  
  compare_table0 <- full_join(compare_table0,shortlist1, by="Biomarkers")
  colnames(compare_table0) <-   c( "n_med", "Biomarkers", "HR_med", "P_med", "P_chk_med", "n_mean", "HR_mean", "P_mean", "P_chk_mean", "n_PH", "HR_PH", "P_PH", "P_chk_PH")
  
  compare_table0 <- full_join(compare_table0,shortlist2, by="Biomarkers")
  
  colnames(compare_table0) <- c( "n_med", "Biomarkers", "HR_med", "P_med", "P_chk_med", "n_mean", "HR_mean", "P_mean", "P_chk_mean", "n_PH", "HR_PH", "P_PH", "P_chk_PH", "n_CM", "HR_CM", "P_CM", "P_chk_CM")
  compare_table0 <- compare_table0 %>%  `rownames<-`( NULL ) %>% column_to_rownames("Biomarkers") %>% mutate_all(~as.numeric(as.character(.)))  %>%  mutate_all(round, digits = 3)
  write.csv(compare_table0, paste0("Cox regression outputs/",name,"_fullMergedTable.csv"), append = TRUE)
  
  ## First comparison, 
  # filter out all significant cox.zph p value 
  # Then show table with either of the three models significant 
  # since this is a large model it is exported as an external table 
  compare_table <- compare_table0 %>% dplyr::filter(`P_mean`<=significance |`P_med`<=significance | `P_CM` <= significance) %>% dplyr::filter(`P_chk_PH`> significance | `P_chk_CM`> significance | `P_chk_med`> significance | `P_chk_mean`> significance )
  
  compare_table %>% select("n_med", "HR_med", "P_med", "n_mean", "HR_mean", "P_mean", "n_PH", "HR_PH", "P_PH", "n_CM", "HR_CM", "P_CM") %>% mutate(
    P_med = ifelse(P_med <= significance,
                   cell_spec(P_med, "html", color = "red", bold = T), 
                   cell_spec(P_med, "html", color = "black", bold = F)), 
    P_mean = ifelse(P_mean <= significance,
                    cell_spec(P_mean, "html", color = "red", bold = T), 
                    cell_spec(P_mean, "html", color = "black", bold = F)), 
    P_PH = ifelse(P_PH <= significance,
                  cell_spec(P_PH, "html", color = "red", bold = T), 
                  cell_spec(P_PH, "html", color = "black", bold = F)), 
    P_CM = ifelse(P_CM <= significance,
                  cell_spec(P_CM, "html", color = "red", bold = T), 
                  cell_spec(P_CM, "html", color = "black", bold = F))) %>%
    kable("html", escape = F, caption = "Significance filtered (OR) condition - All models") %>%
    kable_styling("hover", full_width = F) %>% print()
  
  assign(paste0("CoxmeTableFor_", name), compare_table,  .GlobalEnv)
  
  ## Second comparison,
  # filter out all significant cox.zph p value
  # Then show table with aggregated model as significant

  compare_table2 <- compare_table0 %>% dplyr::filter(`P_mean`<=significance |`P_med`<=significance) %>% dplyr::filter(`P_chk_PH`> significance | `P_chk_CM`> significance | `P_chk_med`> significance | `P_chk_mean`> significance )

  compare_table2%>% select("n_med", "HR_med", "P_med", "n_mean", "HR_mean", "P_mean", "n_PH", "HR_PH", "P_PH", "n_CM", "HR_CM", "P_CM") %>%mutate(
    P_med = ifelse(P_med <= significance,
                   cell_spec(P_med, "html", color = "red", bold = T),
                   cell_spec(P_med, "html", color = "black", bold = F)),
    P_mean = ifelse(P_mean <= significance,
                    cell_spec(P_mean, "html", color = "red", bold = T),
                    cell_spec(P_mean, "html", color = "black", bold = F)),
    P_PH = ifelse(P_PH <= significance,
                  cell_spec(P_PH, "html", color = "red", bold = T),
                  cell_spec(P_PH, "html", color = "black", bold = F)),
    P_CM = ifelse(P_CM <= significance,
                  cell_spec(P_CM, "html", color = "red", bold = T),
                  cell_spec(P_CM, "html", color = "black", bold = F))) %>%
    kable("html", escape = F, caption = "Significance filtered Aggregated models") %>%
    kable_styling("hover", full_width = F)  %>% print()

  ## Third comparison,
  # filter out all significant cox.zph p value
  # Then show table with aggregated model as significant

  compare_table3 <- compare_table0  %>% dplyr::filter(`P_CM`<=significance) %>% dplyr::filter(`P_chk_PH`> significance | `P_chk_CM`> significance | `P_chk_med`> significance | `P_chk_mean`> significance )

  compare_table3 %>% select("n_med", "HR_med", "P_med", "n_mean", "HR_mean", "P_mean", "n_PH", "HR_PH", "P_PH", "n_CM", "HR_CM", "P_CM") %>%mutate(
    P_med = ifelse(P_med <= significance,
                   cell_spec(P_med, "html", color = "red", bold = T),
                   cell_spec(P_med, "html", color = "black", bold = F)),
    P_mean = ifelse(P_mean <= significance,
                    cell_spec(P_mean, "html", color = "red", bold = T),
                    cell_spec(P_mean, "html", color = "black", bold = F)),
    P_PH = ifelse(P_PH <= significance,
                  cell_spec(P_PH, "html", color = "red", bold = T),
                  cell_spec(P_PH, "html", color = "black", bold = F)),
    P_CM = ifelse(P_CM <= significance,
                  cell_spec(P_CM, "html", color = "red", bold = T),
                  cell_spec(P_CM, "html", color = "black", bold = F))) %>%
    kable("html", escape = F, caption = "Significance filtered CoxME") %>%
    kable_styling("hover", full_width = F)  %>% print()

  ## Running multivariate

  if (printforest == TRUE) {
    multivariatetype <- match.arg(multivariatetype)
    multi_table <- data.frame()

    if (multivariatetype == "CoxME") {

      Model_output <- final_uni_coxme_table %>%  dplyr::filter(Pvalue_check <= significance)  %>% select(Biomarkers)

      filtered_table <- final_uni_coxme_table[!(rownames(final_uni_coxme_table) %in% Model_output),]
      filtered_table <- dplyr::filter(filtered_table, Pvalue<=significance)
      if (nrow(filtered_table)==0) {
        print("No significant marker")
        # Here we calculate multivariate coxme model

      }else{
        vl <- filtered_table$Biomarkers
        len <- length(vl)
        subset_markers <- filtered_table$Biomarkers
        subset_df <- cbind(data[,1:s-1], data[subset_markers])
        f1 = as.formula(paste( "Surv(",Time,",", Event, ") ~", paste(unlist(vl), collapse="+"), "+ (1|", PATID,")"))
        md <- coxme(f1, data=subset_df)
        obj <- univ_cox_me(md)
        multi_table <- multi_table %>% rbind(obj) %>% rownames_to_column("Biomarkers")
        colnames(multi_table) <- c("Biomarkers","nvar", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "SE", "Z", "Pvalue","Pvalue_check")
      }
    }else if (multivariatetype == "CoxPH"){
      # Here we calculate multivariate coxPH model
      Model_output <- final_uni_coxph_table %>%  dplyr::filter(Pvalue_check <= significance)  %>% select(Biomarkers)
      filtered_table <- final_uni_coxph_table[!(rownames(final_uni_coxph_table) %in% Model_output),]
      filtered_table <- dplyr::filter(filtered_table, Pvalue<=significance)
      if (nrow(filtered_table)==0) {
        print("No significant marker")

      }else{
        vl <- filtered_table$Biomarkers
        len <- length(vl)
        subset_markers <- filtered_table$Biomarkers
        subset_df <- cbind(data[,1:s-1], data[subset_markers])
        f1 = as.formula(paste("Surv(as.numeric(as.character(subset_df$", Time, ")), as.numeric(as.character(subset_df$", Event, "))) ~", paste(unlist(vl), collapse="+")))
        md <- coxph(f1, data=subset_df)
        obj <- multi_cox_ph(md)
        multi_table <- multi_table %>% rbind(obj)
        multi_table <- multi_table %>% rownames_to_column("Biomarkers")
        colnames(multi_table) <- c("Biomarkers","Vname",  "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
        multi_table <- multi_table %>% select("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
      }
    }else if (multivariatetype == "Aggregated Mean"){
      # Here we calculate multivariate aggreagated model
      Model_output <- final_uni_aggregated_mean_table %>%  dplyr::filter(Pvalue_check <= significance)  %>% select(Biomarkers)
      filtered_table <- final_uni_aggregated_mean_table[!(rownames(final_uni_aggregated_mean_table) %in% Model_output),]
      filtered_table <- dplyr::filter(filtered_table, Pvalue<=significance)

      if (nrow(filtered_table)==0) {
        print("No significant marker")

      }else{
        vl <- filtered_table$Biomarkers
        len <- length(vl)
        subset_markers <- filtered_table$Biomarkers
        subset_df <- cbind(data[,1:s-1], data[subset_markers])
        subset_df <-  subset_df %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, mean)
        f1 = as.formula(paste("Surv(as.numeric(as.character(subset_df$", Time, ")), as.numeric(as.character(subset_df$", Event, "))) ~", paste(unlist(vl), collapse="+")))
        md <- coxph(f1, data=subset_df)
        obj <- multi_cox_ph(md)
        multi_table <- multi_table %>% rbind(obj)
        multi_table <- multi_table %>% rownames_to_column("Biomarkers")
        colnames(multi_table) <- c("Biomarkers","Vname",  "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
        multi_table <- multi_table %>% select("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
      }
    } else if (multivariatetype == "Aggregated Median"){
      # Here we calculate multivariate aggreagated model
      Model_output <- final_uni_aggregated_median_table %>%  dplyr::filter(Pvalue_check <= significance)  %>% select(Biomarkers)
      filtered_table <- final_uni_aggregated_median_table[!(rownames(final_uni_aggregated_median_table) %in% Model_output),]

      filtered_table <- dplyr::filter(filtered_table, Pvalue<=significance)
      if (nrow(filtered_table)==0) {
        print("No significant marker")

      }else{
        vl <- filtered_table$Biomarkers
        len <- length(vl)
        subset_markers <- filtered_table$Biomarkers
        names = vl[j]
        subset_df <- cbind(data[,1:s-1], data[subset_markers])
        subset_df <-  subset_df %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, median)
        f1 = as.formula(paste("Surv(as.numeric(as.character(subset_df$", Time, ")), as.numeric(as.character(subset_df$", Event, "))) ~", paste(unlist(vl), collapse="+")))
        md <- coxph(f1, data=subset_df)
        obj <- multi_cox_ph(md)
        multi_table <- multi_table %>% rbind(obj)
        multi_table <- multi_table %>% rownames_to_column("Biomarkers")
        colnames(multi_table) <- c("Biomarkers","Vname",  "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
        multi_table <- multi_table %>% select("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")

      }
    }else {
      print("No significant marker")
    }

    if (is.null(dim(multi_table))) {
      print("No significant marker")

    } else if (nrow(multi_table)> 0) {
      filtered_final_table2 <- inner_join(filtered_table, multi_table, by = "Biomarkers")

      sfrac <- function(top,bottom,data=NULL) with(data,lapply(paste0("atop(",top,",",bottom,")"),str2expression))
      # Text on plot

      # Plot
      tabletext2<- list(
        c("Biomarkers",filtered_final_table2$Biomarkers),
        c("n",sfrac(filtered_final_table2$n.x,filtered_final_table2$n.y,data=filtered_final_table2)),
        c("HR",sfrac(filtered_final_table2$HR.x,filtered_final_table2$HR.y,data=filtered_final_table2)),
        c(paste0(expression(beta)) ,sfrac(filtered_final_table2$Beta.x,filtered_final_table2$Beta.y, data=filtered_final_table2)),
        c("P",sfrac(filtered_final_table2$Pvalue.x,filtered_final_table2$Pvalue.y,data=filtered_final_table2)))

      pdf(file=paste0(outputfilepath,"CoxMEplot_",name, ".pdf"),width = 18, height = 32 ,onefile = T)
      fp <- forestplot(tabletext2, mean = cbind(c(NA, as.numeric(as.character(filtered_final_table2$HR.x))), c(NA,as.numeric(as.character(filtered_final_table2$HR.y)))),
                       lower = cbind (c(NA,as.numeric(as.character(filtered_final_table2$`95% Lower HR.x`))), c(NA,as.numeric(as.character(filtered_final_table2$`95% Lower HR.y`)))),
                       upper = cbind(c(NA,as.numeric(as.character(filtered_final_table2$`95% Upper HR.x`))), c(NA, as.numeric(as.character(filtered_final_table2$`95% Upper HR.y`)))),
                       clip = c(0.1,5),
                       lineheight = unit(13,"mm"),
                       line.margin = 0.15,
                       xlog = TRUE, xlab = "HR with 95% CI",
                       col = fpColors(box = c("royalblue2", "red3"),
                                      lines = c("royalblue2", "red3")),
                       fn.ci_norm = c(fpDrawCircleCI, fpDrawDiamondCI),
                       #is.summary = c(TRUE,rep(FALSE,5)),
                       graph.pos = 2,
                       boxsize = 0.2,
                       title = paste0(name),
                       xticks = c(0.5, 1, 1.5, 2, 3),
                       legend = c("Univariate", "Multivariate"),
                       vertices = TRUE, txt_gp = fpTxtGp(label=gpar(cex=0.6), ticks = gpar(cex=0.5), xlab = gpar(cex=0.5)), ci.vertices = TRUE, align = c("l","c","c","c", "c","c","c"), graphwidth = unit(4, "cm"), colgap=unit(4,"mm"))
      EnvStats::print(fp)
      dev.off()
      EnvStats::print(fp)
    }
  }
  return(compare_table)
}


KMplotmedian <- function(inputdata, i, aggregationType, aggregate = TRUE, OS ="OS", OS.01="OS.01"){
  
  med <-  median(as.numeric(as.character(inputdata[[i]])), na.rm = TRUE)
  inputdata <- inputdata %>% mutate(Cat = ifelse(inputdata[[i]]<med, "Low", "High"))
  fit <- survfit(Surv(as.numeric(as.character(OS)), as.numeric(as.character(OS.01))) ~ Cat, data = inputdata)
  
  plot <- survminer::ggsurvplot(fit, data = inputdata, pval = TRUE, #conf.int = TRUE,
                                risk.table = TRUE, # Add risk table
                                risk.table.col = "strata", # Change risk table color by groups
                                linetype = "strata", # Change line type by groups
                                surv.median.line = "hv", # Specify median survival
                                ggtheme = theme_bw(), # Change ggplot2 theme
                                title = paste0("Aggregate ", aggregationType, " ", i), 
                                legend.labs =c("High", "Low"), 
                                palette = c("#02818a", "#ae017e"),font.x = 10, font.y = 10,font.legend = 10, font.title=8)
  return(plot)
}

KMCatCurves <- function(data, inputlist, filter=FALSE, filteron, filterCV = FALSE, PATID, CVcutoff = 30, aggregationType = c(NULL, "mean", "median"), SurvCutoff = TRUE, minprop = 0.1, OS ="OS", OS.01="OS.01"){
  
  for (i in inputlist) {
    data[[i]] <- as.numeric(as.character(data[[i]]))
    
    if(filter==TRUE){
      data <- data %>%  dplyr::filter(SegmentLabel == paste0(filteron)) %>% as.data.frame()
      #return(data)
    }else{
      data <- data
      #return(data)
    }
    
    if (filterCV == TRUE) {
      cvlist <- filterCV(data=data, PATID = PATID, names=paste0(i))
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
    }
    
    aggregationType <- match.arg(aggregationType)
    
    if (aggregationType == "median") {
      
      
      datanew <- data%>% select(PATID, OS, OS.01, paste0(i)) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, median) %>% as.data.frame()
      plot1 <- KMplotmedian(inputdata=datanew, i=paste0(i), aggregationType = aggregationType)
      
      
      if (SurvCutoff == TRUE) {
        splots <- list()
        cutpoint <- surv_cutpoint(datanew,time = OS,event = OS.01,variables = paste0(i), minprop = minprop)
        plot(cutpoint) %>% print()
        res.cat <- surv_categorize(cutpoint)
        
        res.cat <- data.frame(res.cat)
        res.cat[[i]] <- as.factor(res.cat[[i]])
        res.cat<-res.cat[complete.cases(res.cat),]
        
        f1 <- as.formula(paste0("Surv(OS, OS.01) ~", i))
        fit <- do.call(survfit, args = list(formula = f1, data = res.cat))
        
        
        plot2 <- survminer::ggsurvplot(fit, data = res.cat, pval = TRUE, conf.int = TRUE,
                                       risk.table = TRUE, # Add risk table
                                       risk.table.col = "strata", # Change risk table color by groups
                                       linetype = "strata", # Change line type by groups
                                       surv.median.line = "hv", # Specify median survival
                                       ggtheme = theme_bw(), # Change ggplot2 theme
                                       title = paste0("Based on ", minprop, " minprop-Survcutoff"),
                                       palette = c("#02818a", "#ae017e"),
                                       legend.labs = c("High", "Low"), font.x = 10, font.y = 10,font.legend = 10, font.title=8)
        
        splots[[1]] <- plot1
        splots[[2]] <- plot2
        
        s <- arrange_ggsurvplots(splots,print = TRUE,title = NA,ncol = 2,nrow = 1,surv.plot.height = NULL, risk.table.height = NULL,ncensor.plot.height = NULL)
        
        
      }else {
        print(plot1)
      }
      
    } else if (aggregationType == "mean") {
      
      datanew <- data %>% select(PATID, OS, OS.01, paste0(i)) %>%  dplyr::group_by(PATID) %>% summarise_if(is.numeric, mean) %>% as.data.frame()
      plot1 <-KMplotmedian(inputdata=datanew, i=paste0(i), aggregationType = aggregationType)
      
      if (SurvCutoff == TRUE) {
        splots <- list()
        cutpoint <- surv_cutpoint(datanew,time = OS,event = OS.01,variables = paste0(i), minprop = minprop)
        plot(cutpoint) %>% print()
        res.cat <- surv_categorize(cutpoint)
        res.cat <- data.frame(res.cat)
        res.cat[[i]] <- as.factor(res.cat[[i]])
        res.cat<-res.cat[complete.cases(res.cat),]
        
        f1 <- as.formula(paste0("Surv(OS, OS.01) ~", i))
        fit <- do.call(survfit, args = list(formula = f1, data = res.cat))
        
        plot2 <- survminer::ggsurvplot(fit, data = res.cat, pval = TRUE, conf.int = TRUE,
                                       risk.table = TRUE, # Add risk table
                                       risk.table.col = "strata", # Change risk table color by groups
                                       linetype = "strata", # Change line type by groups
                                       surv.median.line = "hv", # Specify median survival
                                       ggtheme = theme_bw(), # Change ggplot2 theme
                                       title = paste0("Based on ", minprop, " minprop-Survcutoff"),
                                       palette = c("#02818a", "#ae017e"),
                                       legend.labs = c("High", "Low"),  font.x = 10, font.y = 10,font.legend = 10, font.title=8)
        
        splots[[1]] <- plot1
        splots[[2]] <- plot2
        
        s <- arrange_ggsurvplots(splots,print = TRUE,title = NA,ncol = 2,nrow = 1,surv.plot.height = NULL, risk.table.height = NULL,ncensor.plot.height = NULL)
        
        
        
      }else {
        print(plot1)
      }
      
      
    } else if (is.null(aggregationType)) {
      plot1 <-KMplotmedian(inputdata=data, i=paste0(i), aggregationType = "No aggregation")
      
      if (SurvCutoff == TRUE) {
        splots <- list()
        cutpoint <- surv_cutpoint(data,time = OS,event = OS.01,variables = paste0(i), minprop = minprop)
        plot(cutpoint) %>% print()
        res.cat <- surv_categorize(cutpoint)
        #fit <- survfit(Surv(as.numeric(as.character(OS)), as.numeric(as.character(OS.01))) ~ paste0(i), data = res.cat)
        res.cat <- data.frame(res.cat)
        res.cat[[i]] <- as.factor(res.cat[[i]])
        res.cat<-res.cat[complete.cases(res.cat),]
        
        f1 <- as.formula(paste0("Surv(OS, OS.01) ~", i))
        fit <- do.call(survfit, args = list(formula = f1, data = res.cat))
        
        plot2 <- survminer::ggsurvplot(fit, data = res.cat, pval = TRUE, conf.int = TRUE,
                                       risk.table = TRUE, # Add risk table
                                       risk.table.col = "strata", # Change risk table color by groups
                                       linetype = "strata", # Change line type by groups
                                       surv.median.line = "hv", # Specify median survival
                                       ggtheme = theme_bw(), # Change ggplot2 theme
                                       title = paste0("Based on ", minprop, " minprop-Survcutoff"),
                                       palette = c("#02818a", "#ae017e"),
                                       legend.labs = c("High", "Low"), font.x = 10, font.y = 10,font.legend = 10, font.title=8)
        
        splots[[1]] <- plot1
        splots[[2]] <- plot2
        
        s <- arrange_ggsurvplots(splots,print = TRUE,title = NA,ncol = 2,nrow = 1,surv.plot.height = NULL, risk.table.height = NULL,ncensor.plot.height = NULL)
        
        
        
      }else {
        print(plot1)
      }
    }
  }
}

coxphfunction <- function(data, name, s, e, filterCV=TRUE, CVcutoff, Time, Event){
  vl <- names(data[,s:e])
  len <- length(vl)
  final_uni_coxph_table <- data.frame()
  
  # This chunk calculates coxph
  for (j in 1:len){
    if (filterCV == TRUE) {
      #print(vl[j])
      cvlist <- filterCV(data=data, PATID = "PATID", names=vl[j])
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data2 <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
      
      output <- RunCoxPH(data = data2, names=vl[j], Time = Time, Event=Event)
      final_uni_coxph_table <- final_uni_coxph_table %>% rbind(output)
      
    } else{
      output <- RunCoxPH(data = data, names=vl[j], Time = Time, Event= Event)
      final_uni_coxph_table <- final_uni_coxph_table %>% rbind(output)
    }
  }
  
  colnames(final_uni_coxph_table) <- c("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
  write.csv(final_uni_coxph_table, paste0("Cox regression outputs/",name,"_coxPH.csv"), append = TRUE)
  return(final_uni_coxph_table)
}

coxmefunction <- function(data, name, s, e, filterCV=TRUE, CVcutoff, Time, Event){
  vl <- names(data[,s:e])
  len <- length(vl)
  final_uni_coxme_table <- data.frame()
  for (j in 1:len){
    names = vl[j]
    if (filterCV == TRUE) {
      
      cvlist <- filterCV(data=data,  PATID = "PATID", names=vl[j])
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data2 <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
      f1 = as.formula(paste("Surv(", Time, ",", Event,") ~", paste(vl[j], "+ (1|PATID)")))
      skip_to_next <- FALSE
      tryCatch(coxme(f1, data=data2), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { next }
      md <- coxme(f1, data2)
      obj <- univ_cox_me(md)
      final_uni_coxme_table <- final_uni_coxme_table %>% rbind(obj)
      
    } else {
      
      f1 = as.formula(paste("Surv(", Time, ",", Event,") ~", paste(vl[j], "+ (1|PATID)")))
      skip_to_next <- FALSE
      tryCatch(coxme(f1, data=data), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { next }
      md <- coxme(f1, data)
      obj <- univ_cox_me(md)
      final_uni_coxme_table <- final_uni_coxme_table %>% rbind(obj)
    }
  }
  colnames(final_uni_coxme_table) <- c("nvar", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "SE", "Z", "Pvalue","Pvalue_check")
  final_uni_coxme_table <- final_uni_coxme_table %>%rownames_to_column("Biomarkers") 
  write.csv(final_uni_coxme_table, paste0("Cox regression outputs/",name,"_coxME.csv"), append = TRUE)
  return(final_uni_coxme_table)
}

Mean_A_function <- function(data, name, s, e, filterCV=TRUE, CVcutoff, Time, Event, uniquefactortype = NULL){
  
  vl <- names(data[,s:e])
  len <- length(vl)
  final_uni_aggregated_mean_table <- data.frame()
  
  for (j in 1:len){
    names = vl[j]
    
    if (filterCV == TRUE) {
      
      cvlist <- filterCV(data=data,  PATID = "PATID", names=names)
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data2 <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
      
      if (is.null(uniquefactortype)) {
        shortdf <-  data2 %>% dplyr::select(PATID,Time, Event, vl[j]) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, mean)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_mean_table <- final_uni_aggregated_mean_table %>% rbind(output)
        # print(vl[j])
        
      } else {
        shortdf <-  data2 %>% dplyr::select(PATID, paste0(uniquefactortype),Time, Event,  vl[j]) %>% dplyr::group_by(PATID, paste0(uniquefactortype)) %>% summarise_if(is.numeric, mean)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_mean_table <- final_uni_aggregated_mean_table %>% rbind(output)
        
      }
    } else {
      if (is.null(uniquefactortype)) {
        shortdf <-  data %>% dplyr::select(PATID, Time, Event, vl[j]) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, mean)
        
        output <- RunCoxPH(data = shortdf, names=vl[j], Time = "OS", Event="OS.01")
        final_uni_aggregated_mean_table <- final_uni_aggregated_mean_table %>% rbind(output)
        
      } else {
        shortdf <-  data %>% dplyr::select(PATID,paste0(uniquefactortype), Time, Event, vl[j]) %>% dplyr::group_by(PATID, paste0(uniquefactortype))  %>% summarise_if(is.numeric, mean)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_mean_table <- final_uni_aggregated_mean_table %>% rbind(output)
        
      }
    }
  }
  
  colnames(final_uni_aggregated_mean_table) <- c("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
  write.csv(final_uni_aggregated_mean_table, paste0("Cox regression outputs/",name,"_A-Mean.csv"), append = TRUE)
  return(final_uni_aggregated_mean_table)
}

Median_A_function <- function(data, name, s, e, filterCV=TRUE, CVcutoff, Time, Event, uniquefactortype = NULL){
  
  vl <- names(data[,s:e])
  len <- length(vl)
  final_uni_aggregated_median_table <- data.frame()
  
  for (j in 1:len){
    names = vl[j]
    
    if (filterCV == TRUE) {
      
      cvlist <- filterCV(data=data,  PATID = "PATID", names=names)
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data2 <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
      
      if (is.null(uniquefactortype)) {
        shortdf <-  data2 %>% dplyr::select(PATID,Time, Event, vl[j]) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, median)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_median_table <- final_uni_aggregated_median_table %>% rbind(output)
        # print(vl[j])
        
      } else {
        shortdf <-  data2 %>% dplyr::select(PATID, paste0(uniquefactortype),Time, Event,  vl[j]) %>% dplyr::group_by(PATID, paste0(uniquefactortype)) %>% summarise_if(is.numeric, median)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_median_table <- final_uni_aggregated_median_table %>% rbind(output)
        
      }
    } else {
      if (is.null(uniquefactortype)) {
        shortdf <-  data %>% dplyr::select(PATID, Time, Event, vl[j]) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, median)
        
        output <- RunCoxPH(data = shortdf, names=vl[j], Time = "OS", Event="OS.01")
        final_uni_aggregated_median_table <- final_uni_aggregated_median_table %>% rbind(output)
        
      } else {
        shortdf <-  data %>% dplyr::select(PATID,paste0(uniquefactortype), Time, Event, vl[j]) %>% dplyr::group_by(PATID, paste0(uniquefactortype))  %>% summarise_if(is.numeric, median)
        
        output <- RunCoxPH(data = shortdf, names=names, Time = "OS", Event="OS.01")
        final_uni_aggregated_median_table <- final_uni_aggregated_median_table %>% rbind(output)
        
      }
    }
  }
  
  colnames(final_uni_aggregated_median_table) <- c("Biomarkers", "n", "Beta","HR", "95% Lower HR", "95% Upper HR", "Pvalue","Pvalue_check")
  write.csv(final_uni_aggregated_median_table, paste0("Cox regression outputs/",name,"_A-median.csv"), append = TRUE)
  return(final_uni_aggregated_median_table)
}

logrankmedianfilter <- function(inputdata, i){
  
  med <-  median(as.numeric(as.character(inputdata[[i]])), na.rm = TRUE)
  inputdata <- inputdata %>% mutate(Cat = ifelse(inputdata[[i]]<med, "Low", "High"))
  fit <- survdiff(Surv(as.numeric(as.character(OS)), as.numeric(as.character(OS.01))) ~ Cat, data = inputdata)
  p.val <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
  
  return(p.val)
}


logrankFilter <- function(data, inputlist, filterCV = FALSE, PATID, CVcutoff = 30, aggregationType = c(NULL, "mean", "median"), SurvCutoff = TRUE, minprop = 0.1){
  splots <- data.frame()
  for (i in inputlist) {
    data[[i]] <- as.numeric(as.character(data[[i]]))
    
    
    if (filterCV == TRUE) {
      cvlist <- filterCV(data=data, PATID = PATID, names=paste0(i))
      cvlist <- cvlist %>%  dplyr::filter(CV <= CVcutoff)
      
      data <- subset(data, data$PATID %in% cvlist$PATID) %>% as.data.frame()
    }
    
    aggregationType <- match.arg(aggregationType)
    
    if (aggregationType == "median") {
      
      datanew <- data%>% select(PATID, OS, OS.01, paste0(i)) %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, median) %>% as.data.frame()
      Median_cutoff <- logrankmedianfilter(inputdata=datanew, i=paste0(i))
      
      if (SurvCutoff == TRUE) {
        
        cutpoint <- surv_cutpoint(datanew,time = "OS",event = "OS.01",variables = paste0(i), minprop = minprop)
        res.cat <- surv_categorize(cutpoint)
        res.cat[[i]] <- as.factor(res.cat[[i]])
        res.cat<-res.cat[complete.cases(res.cat),]
        
        f1 <- as.formula(paste0("Surv(OS, OS.01) ~", i))
        fit <- do.call(survdiff, args = list(formula = f1, data = res.cat))
        Survcutoff <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
        
        collated <-  cbind(paste0(i), Median_cutoff, Survcutoff)
        splots <- splots %>% rbind(collated)  
        
        
      }else {
        collated <-  cbind(paste0(i), Median_cutoff)
        splots <- splots %>% rbind(collated)
      }
      
    } else if (aggregationType == "mean") {
      
      datanew <- data %>% select(PATID, OS, OS.01, paste0(i)) %>%  dplyr::group_by(PATID) %>% summarise_if(is.numeric, mean) %>% as.data.frame()
      Median_cutoff <- logrankmedianfilter(inputdata=datanew, i=paste0(i))
      
      if (SurvCutoff == TRUE) {
        splots <- list()
        cutpoint <- surv_cutpoint(datanew,time = "OS",event = "OS.01",variables = paste0(i), minprop = minprop)
        
        res.cat <- surv_categorize(cutpoint)
        res.cat[[i]] <- as.factor(res.cat[[i]])
        res.cat<-res.cat[complete.cases(res.cat),]
        
        f1 <- as.formula(paste0("Surv(OS, OS.01) ~", i))
        fit <- do.call(survdiff, args = list(formula = f1, data = res.cat))
        Survcutoff <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
          
        collated <-  cbind(paste0(i), Median_cutoff, Survcutoff)
        splots <- splots %>% rbind(collated)  
        
      }else {
        collated <-  cbind(paste0(i), Median_cutoff)
        splots <- splots %>% rbind(collated)
      }
      
      
    } else if (is.null(aggregationType)) {
      Median_cutoff <- logrankmedianfilter(inputdata=datanew, i=paste0(i))
      
      if (SurvCutoff == TRUE) {
        splots <- list()
        cutpoint <- surv_cutpoint(data,time = "OS",event = "OS.01",variables = paste0(i), minprop = minprop)
        
        res.cat <- surv_categorize(cutpoint)
        res.cat[[i]] <- as.factor(res.cat[[i]])
        res.cat<-res.cat[complete.cases(res.cat),]
        
        f1 <- as.formula(paste0("Surv(OS, OS.01) ~", i))
        fit <- do.call(survdiff, args = list(formula = f1, data = res.cat))
        Survcutoff <- 1 - pchisq(fit$chisq, length(fit$n) - 1)
        collated <-  cbind(paste0(i), Median_cutoff, Survcutoff)
        splots <- splots %>% rbind(collated)  
        
      }else {
        collated <-  cbind(paste0(i), Median_cutoff)
        splots <- splots %>% rbind(collated)
      }
    }
  }
  return(splots)
}
