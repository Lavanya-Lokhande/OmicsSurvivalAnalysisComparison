KMCatCurves <- function(data, inputlist, filter=TRUE, filteron){
  for (i in inputlist) {
    print(i)
    data[[i]] <- as.numeric(as.character(data[[i]]))
    head(data[[i]])
    
    if(filter==TRUE){
      data <- data %>%  dplyr::filter(SegmentLabel == paste0(filteron)) %>% as.data.frame()
      #return(data)
    }else{
      data <- data
      #return(data)
    }
    
    med <-  median(as.numeric(as.character(data[[i]])), na.rm = TRUE)
    
    data <- data %>% mutate(Cat = ifelse(data[[i]]<med, "Low", "High"))
    
    
    fit <- survfit(Surv(as.numeric(as.character(OS)), as.numeric(as.character(OS.01))) ~ Cat, data = data)
    
    plot <- survminer::ggsurvplot(fit, data = data, pval = TRUE, conf.int = TRUE,
                                  risk.table = TRUE, # Add risk table
                                  risk.table.col = "strata", # Change risk table color by groups
                                  linetype = "strata", # Change line type by groups
                                  surv.median.line = "hv", # Specify median survival
                                  ggtheme = theme_bw(), # Change ggplot2 theme
                                  title = paste0("All_", i))
    EnvStats::print(plot)
    
    datanew <- data %>% dplyr::group_by(PATID) %>% summarise_if(is.numeric, mean) %>% as.data.frame()
    med <-  median(as.numeric(as.character(datanew[[i]])), na.rm = TRUE)
    datanew <- datanew %>% mutate(Cat2 = ifelse(datanew[[i]]<med, "Low", "High"))
    
    fit <- survfit(Surv(as.numeric(as.character(OS)), as.numeric(as.character(OS.01))) ~ Cat2, data = datanew)
    
    plot2 <- survminer::ggsurvplot(fit, data = datanew, pval = TRUE, conf.int = TRUE,
                                   risk.table = TRUE, # Add risk table
                                   risk.table.col = "strata", # Change risk table color by groups
                                   linetype = "strata", # Change line type by groups
                                   surv.median.line = "hv", # Specify median survival
                                   ggtheme = theme_bw(), # Change ggplot2 theme
                                   title = paste0("PATID_averaged_", i))
    
    EnvStats::print(plot2)
  }
  
}