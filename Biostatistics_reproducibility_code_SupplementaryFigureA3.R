library(forestplot)
library(grid)
library(metaRMST)
library(survRM2)

# first obtain effect estimates for each trial at each time horizon


# determine RMSTD at each time point:
rmst_calculations <- function(time_horizons, data, trialnames){

  nrows <- max(data$trialID)*length(time_horizons)
  result <- data.frame(study=numeric(nrows), time_horizon=numeric(nrows), Estimate=numeric(nrows), SE=numeric(nrows), lower=numeric(nrows), upper=numeric(nrows))

  for (j in 1:max(data$trialID)){
    index <- -max(data$trialID)+j
    dat <- AorticStenosisTrials[which(AorticStenosisTrials$trialID==j),]

    FU <- min(max(dat[which(dat$Arm==0),]$Time), max(dat[which(dat$Arm==1),]$Time)) # minimum of max observed followup times across groups

    for(tau in time_horizons){
      index <- index+max(data$trialID)

      if(FU>=tau){obj<-rmst2(dat$Time, dat$Event, dat$Arm, tau=tau)}else{obj <- NA}
      result$time_horizon[index] <- tau
      result$study[index] <- trialnames[j]
      if(FU>=tau){result$Estimate[index] <- obj$unadjusted.result[1]}else{result$Estimate[index] <- NA}
      if(FU>=tau){result$SE[index]  <-  (obj$unadjusted.result[1,2]-obj$unadjusted.result[1])/-qnorm(0.975)}else{result$SE[index]  <- NA}
      if(FU>=tau){result$lower[index] <-  obj$unadjusted.result[1,2]}else{result$lower[index] <- NA}
      if(FU>=tau){result$upper[index] <-  obj$unadjusted.result[1,3]}else{result$upper[index] <- NA}

    }
  }
  return(result)
}

rmst_res <- rmst_calculations(c(12,24, 36), AorticStenosisTrials, trialnames=c("Notion", "Partner", "Surtavi", "Partner2", "Corevalve"))
rmst_res <- rmst_res[order(rmst_res$time_horizon),]

# then obtain meta-analytic combined effect estimates at each time horizon
mvma_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="mvma")
boot_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="mvma_boot", nboot=500)
uni_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="uni")
uni_flex_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="uni_flex")
results <- rbind(mvma_res$REresult, boot_res$REresult, uni_flex_res$result[,1:9], uni_res$result[,1:9])
results_reduced <- results[,c(1,2,4:7)]
colnames(results_reduced)[1] <- "study"

# combine the results and arrange 
forest_res <- rbind(results_reduced, rmst_res)
forest_res <- forest_res[order(forest_res$time_horizon),]

forest_res$index <- ifelse(forest_res$study=="Partner",1,
	ifelse(forest_res$study=="Corevalve",2,
	ifelse(forest_res$study=="Notion",3,
	ifelse(forest_res$study=="Surtavi",4,
	ifelse(forest_res$study=="Partner2",5,
	ifelse(forest_res$study=="univariate",6,
	ifelse(forest_res$study=="Univariate with model estimates",7,
	ifelse(forest_res$study=="Random Effect MVMA",8,
	ifelse(forest_res$study=="Random Effect MVMA boot",9,999)))))))))

forest_res<- forest_res[order(forest_res$time_horizon, forest_res$index),]

# prepare plot space:
data12 <- forest_res[which(forest_res$time_horizon==12),]
np12 <- paste(format(round(data12$Estimate, digits=2),nsmall=1)," (",format(round(data12$lower,digits=2),nsmall=1),", ",format(round(data12$upper,digits=2),nsmall=1),")",sep="")
tabletext <- cbind(c("Study", "Partner", "Corevalve", "Notion", "Surtavi", "Partner2", "Univariate, available data", "Univariate, flex. model estimates", "Multivariate, analytic covariance", "Multivariate, bootstrap covariance"))

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
pushViewport(viewport(layout.pos.col = 1))

# invisible plot
forestplot(labeltext=tabletext, graph.pos=1, 
           mean=rep(NA,10), 
           lower=rep(NA,10), upper=rep(NA,10),
           #title="",
           xlab=" \n ", is.summary=c(TRUE,  rep(FALSE,5), rep(TRUE, 4)),
           txt_gp=fpTxtGp(label=gpar(cex=1.2, fontface="plain"),
                              ticks=gpar(cex=1.1, col="white"),
                              xlab=gpar(cex = 1.2),
                              title=gpar(cex = 1.2, fontface="plain")),
  		#hrzl_lines=list("2" = gpar(lwd=1, col="black")),
           col=fpColors(box=NA, lines=NA, zero=NA),
           zero=5, xticks=c(0,1), cex=0.9, lineheight = "auto", boxsize=0.1, colgap=unit(1,"mm"),
	   #fn.ci_norm=matrix(c("fpDrawCircleCI"), nrow=17,ncol=1, byrow=T),
           lwd.ci=1, ci.vertices=FALSE, ci.vertices.height = 0, new_page=FALSE)

popViewport()

grid.rect(x = unit(0.01, "npc"), y = unit(0.1, "npc"),
          width = unit(.04, "npc"), height = unit(.1, "npc"),
          just = "centre", hjust = NULL, vjust = NULL,
          default.units = "npc", name = NULL,
          gp=gpar(col="white", fill="white"), draw = TRUE, vp = NULL)

# add 12 month results:
# 12 month plot 
tabletext12 <- cbind(c(" ", np12))

pushViewport(viewport(layout.pos.col = 2))

forestplot(labeltext=tabletext12, graph.pos=1, 
           mean=c(NA,data12$Estimate), 
           lower=c(NA,data12$lower), upper=c(NA,data12$upper),
           #title="Difference in Mean Likert Score", graphwidth=unit(18, "mm"),
           xlab=" \n ", is.summary=c(TRUE,  rep(FALSE,5), rep(TRUE, 4)),
           txt_gp=fpTxtGp(label=gpar(cex=1, fontface="plain"),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2, fontface="plain")),
           #hrzl_lines=list("2" = gpar(lwd=1, col="black"), "7"=gpar(lty=2)),
           col=fpColors(box="black", lines="black", zero="black", hrz_lines="#444444"),
           zero=0, xticks=c(-2, -1,0,1,2,3), cex=0.9, lineheight = "auto", boxsize=0.1, colgap=unit(12,"mm"),
           fn.ci_norm=matrix(c("fpDrawCircleCI"), nrow=10,ncol=1, byrow=T),
           #fn.ci_sum=matrix(c("fpDrawSummaryCI"), nrow=1, ncol=1, byrow=T),
           lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.1, new_page=FALSE)

grid.text('Time Horizon 12 months', 0.29, .92, gp = gpar(fontsize=12, font = 1))
grid.text('RMSTD', 0.22, 0.05, gp = gpar(fontsize=12, font = 1))


popViewport()

# 24 month plot
# 24 month time horizon forest plot: 
pushViewport(viewport(layout.pos.col = 3))

data24 <- forest_res[which(forest_res$time_horizon==24),]
np24 <- paste(format(round(data24$Estimate, digits=2),nsmall=1)," (",format(round(data24$lower,digits=2),nsmall=1),", ",format(round(data24$upper,digits=2),nsmall=1),")",sep="")
#np24[3:5] <- ""
tabletext24 <- cbind(c("", np24))

forestplot(labeltext=tabletext24, graph.pos=1, 
           mean=c(NA,data24$Estimate), 
           lower=c(NA,data24$lower), upper=c(NA,data24$upper),
           #title="Difference in Mean Likert Score",
           xlab=" \n ", is.summary=c(TRUE,  rep(FALSE,5), rep(TRUE, 4)),
           txt_gp=fpTxtGp(label=gpar(cex=1, fontface="plain"),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2, fontface="plain")),
           #hrzl_lines=list("2" = gpar(lwd=1, col="black"), "7"=gpar(lty=2)),
           col=fpColors(box="black", lines="black", zero="black", hrz_lines="#444444"),
           zero=0, xticks=c(-2, -1,0,1,2,3), cex=0.9, lineheight = "auto", boxsize=0.1, colgap=unit(6,"mm"),
           fn.ci_norm=matrix(c("fpDrawCircleCI"), nrow=10,ncol=1, byrow=T),
           #fn.ci_sum=matrix(c("fpDrawSummaryCI"), nrow=1, ncol=1, byrow=T),
           lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.1, new_page=FALSE)

grid.text('Time Horizon 24 months', 0.29, .92, gp = gpar(fontsize=12, font = 1))
grid.text('RMSTD', 0.22, 0.05, gp = gpar(fontsize=12, font = 1))

popViewport()

#36 month plot
pushViewport(viewport(layout.pos.col = 4))

data36 <- forest_res[which(forest_res$time_horizon==36),]
np36 <- paste(format(round(data36$Estimate, digits=2),nsmall=1)," (",format(round(data36$lower,digits=2),nsmall=1),", ",format(round(data36$upper,digits=2),nsmall=1),")",sep="")
np36[3:5] <- ""
tabletext36 <- cbind(c("", np36))
tabletext36[8,1] <- "0.54 (-0.05, 1.12)"

forestplot(labeltext=tabletext36, graph.pos=1, 
           mean=c(NA,data36$Estimate), 
           lower=c(NA,data36$lower), upper=c(NA,data36$upper),
           #title="Difference in Mean Likert Score",
           xlab=" \n ", is.summary=c(TRUE,  rep(FALSE,5), rep(TRUE, 4)),
           txt_gp=fpTxtGp(label=gpar(cex=1, fontface="plain"),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2, fontface="plain")),
           #hrzl_lines=list("2" = gpar(lwd=1, col="black"), "7"=gpar(lty=2)),
           col=fpColors(box="black", lines="black", zero="black", hrz_lines="#444444"),
           zero=0, xticks=c(-2, -1,0,1,2,3), cex=0.9, lineheight = "auto", boxsize=0.1, colgap=unit(6,"mm"),
           fn.ci_norm=matrix(c("fpDrawCircleCI"), nrow=10,ncol=1, byrow=T),
           #fn.ci_sum=matrix(c("fpDrawSummaryCI"), nrow=1, ncol=1, byrow=T),
           lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.1, new_page=FALSE)

grid.text('Time Horizon 36 months', 0.29, .92, gp = gpar(fontsize=12, font = 1))
grid.text('RMSTD', 0.22, 0.05, gp = gpar(fontsize=12, font = 1))

popViewport()

grid.lines(x = unit(c(0.02,.995), "npc"),
          y = unit(c(0.89, 0.89), "npc"),
          default.units = "npc",
          arrow = NULL, name = NULL,
          gp=gpar(col="black"), draw = TRUE, vp = NULL)

grid.lines(x = unit(c(0.02,.995), "npc"),
          y = unit(c(0.45, 0.45), "npc"),
          default.units = "npc",
          arrow = NULL, name = NULL,
          gp=gpar(col="black", lty=2), draw = TRUE, vp = NULL)













