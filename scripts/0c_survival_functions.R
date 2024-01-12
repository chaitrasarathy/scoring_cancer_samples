survivalAnalysis <- function(survTab, sigScoreTab, type, method="n", cutType="n", cutVal=1, plotTitle){
  surv.sig.merge <- survTab %>%
    inner_join(sigScoreTab, by=c("sample"))

  #Optimizing cut off
  if (type == "diapause"){
    if (cutType == "manual"){
      cutoff <- cutVal
      res.cutp <- 0
      # tryCatch(
      #   expr = {
      #     ifelse(min(surv.sig.merge$dia.score) < cutoff,
      #            warning("Proceeding to survival analysis... "),
      #            warning("Min value greater than cut-off"))
      #     },
      #   error = {
      #     message('Caught an error!')
      #     print(e)
      #   },
      #   warning = {
      #   stop(call. = F)
      #   message("Stopping execution.")
      #   }
      # )

    }else{
      res.cutp <- surv_cutpoint(surv.sig.merge, time = "over_surv", event = "vital_status", variables = c("dia.score"))
      print(res.cutp)                                                                  #Check best cut-off
      cutoff <- res.cutp$cutpoint$cutpoint
    }

    val <- min(surv.sig.merge$dia.score) > cutoff
    # if(val==T){
    #   stop("Exiting", silent = T, call. = F, domain = NA)
    # }

    if(is.na(cutoff)){
      print("Median cutoff")
      cutoff <- median(surv.sig.merge$dia.score)
    }

    # if (val){
    surv.sig.merge$select <- ""                                                       #Adding a column for high/low groups, based on previous cutoff
    surv.sig.merge[(surv.sig.merge$dia.score > cutoff), "select"] <- "2.High"
    surv.sig.merge[(surv.sig.merge$dia.score <= cutoff ), "select"] <- "1.Low"

    #survival analysis
    surv.mod <- survfit(Surv(over_surv, vital_status) ~ select, data = surv.sig.merge)
    surv.plot <- ggsurvplot(surv.mod, data = surv.sig.merge,
                            title = plotTitle,
                            palette =  c("#FF9E29", "#86AA00"),
                            # palette = c("#E7B800", "#2E9FDF"),
                            # palette = c("black", "red"),
                            pval = TRUE, conf.int = F,
                            risk.table = TRUE, xlab = "Days", ylab="Probability")

    cox.mod <- coxph(Surv(over_surv, vital_status) ~ select, data = surv.sig.merge)
    # summary(cox.mod)
    # }else{
    #   message("Bye.")
    #   stop(call. = F, domain = NA)
    # }
  }

  if (type == "autophagy"){
    if(method == "n"){
      res.cutp <- surv_cutpoint(surv.sig.merge, time = "over_surv", event = "vital_status", variables = c("ap.score"))

      print(res.cutp)
      cutoff <- res.cutp$cutpoint$cutpoint

      # wang 21 paper - cutoff is different
    }
    if(method == "w"){
      # res.cutp <- surv_cutpoint(surv.sig.merge, time = "over_surv", event = "vital_status", variables = c("ap.score"))

      # print(res.cutp)                                                                  #Check best cut-off
      cutoff <- quantile(survAn.aut$surv.sig.merge$ap.score)[4]

    }
    if (cutType == "manual"){
      cutoff <- cutVal
      res.cutp <- 0
    }

    if(is.na(cutoff)){
      print("Median cutoff")
      cutoff <- median(surv.sig.merge$ap.score)
    }

    surv.sig.merge$select <- 0
    surv.sig.merge[(surv.sig.merge$ap.score > cutoff), "select"] <- "2.High"
    surv.sig.merge[(surv.sig.merge$ap.score <= cutoff ), "select"] <- "1.Low"

    #survival analysis
    surv.mod <- survfit(Surv(over_surv, vital_status) ~ select, data = surv.sig.merge)
    surv.plot <- ggsurvplot(surv.mod, data = surv.sig.merge,
                            title = plotTitle,
                            palette =  c("#FF9E29", "#86AA00"),
                            # palette = c("#E7B800", "#2E9FDF"),
                            # palette = c("black", "red"),
                            pval = TRUE, conf.int = F,
                            risk.table = TRUE, xlab = "Days", ylab="Probability")

    cox.mod <- coxph(Surv(over_surv, vital_status) ~ select, data = surv.sig.merge)

  }

  # print(cutoff)

  return(list("surv.sig.merge" = surv.sig.merge, "res.cutp" = ifelse(any(grep(x = ls(), pattern = "res.cutp")), res.cutp, cutoff), "cutoff"=cutoff, "surv.mod" = surv.mod, "surv.plot" = surv.plot, "cox.mod" = cox.mod))
  # return(list("surv.sig.merge" = surv.sig.merge, "res.cutp" = res.cutp, "surv.mod" = surv.mod))
}
