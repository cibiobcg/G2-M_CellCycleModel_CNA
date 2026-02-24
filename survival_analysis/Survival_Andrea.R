library(ggplot2)
library(survival)
library(survminer)
library(data.table)
library(RColorBrewer)
library(rstudioapi)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

survival <- fread("../Files_generated/TCGATimePoints.tsv",data.table=F)

sims <- fread("../Files_generated/final_patient.csv",data.table=F)
sims$Patient_ID <- sapply(sims$Patient,function(x) paste(strsplit(x,"-")[[1]][1:3],collapse="-"))

survival <- survival[which(survival$type=="BRCA"),]
survival$sim_time <- NA
survival$sim_perc <- NA
survival$subtype <- NA

## subtype

subtypes <- fread("../Files_generated/Cancer_subtypes_list.tsv",data.table=F)

##########

for(i in 1:nrow(survival)){
  t <- which(sims$Patient_ID==survival$bcr_patient_barcode[i])
  if(length(t)>0)
  {
    survival$sim_time[i] <- sims$Mean[t] # for every patient associo il tempo medio a cui avviene lo switch
    survival$sim_perc[i] <- sims$Percentage_switch[t] # for every patients i associate its switch percentage: quante volte su cento avviene lo switch
  }
  t <- which(subtypes$COMMON==survival$bcr_patient_barcode[i])
  if(length(t)>0)
  {
    survival$subtype[i] <- subtypes$Subtype[t][1] # associo anche il subtype
  }
}

survival <- survival[which(!is.na(survival$sim_time)),]


# survival$switch_time <- 0
# survival <- survival[which(!is.na(survival$sim_time)),]
# survival$switch_time[which(survival$sim_time>3290)] = 1 # switch sopra 3000
# table(survival$switch_time) # mette i pazienti dopra 3000 switch time ad uno e sotto a 0 
# 
# survival$switch_perc <- 0
# survival <- survival[which(!is.na(survival$sim_time)),]
# survival$switch_perc[which(survival$sim_perc<55)] = 1
# table(survival$switch_perc) # mette i pazienti sotto 50 switch ad uno e sotto a 0 

hist(survival$sim_time,breaks=20,main="Switch time distribution",xlab="Patient mean switch time")
hist(survival$sim_perc[which(survival$sim_perc<100)],breaks=10,main="Switch percentage distribution\n(only for patients with percentages<100)",xlab="Patient mean switch percentage")
hist(survival$age_at_initial_pathologic_diagnosis,breaks=20,main="Age distribution",xlab="Age")

summary(survival$sim_time)


d<-density(na.omit(survival$sim_time))
d$x[which.min(abs(diff(d$y)))]
v<-optimize(approxfun(d$x, d$y), interval = c(50,9000))$minimum
# we show this by doing:
hist(survival$sim_time,breaks=20,main="Switch time distribution",xlab="Patient mean switch time", prob=T)
#lines(d, col='red', lty=2)
abline(v=3290, col='blue') # trial

km<-kmeans(survival$sim_time, centers = 2)
survival$clust<-as.factor(km$cluster)
ggplot(survival, aes(x=sim_time, fill=clust))+geom_histogram(alpha=0.5, color='grey50', aes(y = after_stat(count / nrow(survival))))+labs(x='Patients mean switch time', y='Percentage of patients', title='Switch time distribution') + scale_fill_manual(name='',labels=c('Delayed','Fast'), values=c('#eeb78c','#98d3c1'))



prova<-survival[order(survival$sim_time),,drop=F]

## For the switch percentage
# km_perc<-kmeans(survival$sim_perc[which(survival$sim_perc<100)], centers = 2)
km_perc<-kmeans(survival$sim_perc, centers = 2)

# prova<-survival[which(survival$sim_perc<100),]
prova <-survival
prova$clust2<-as.factor(km_perc$cluster)
# library(ggplot2)
ggplot(prova, aes(x=sim_perc, fill=clust2))+geom_histogram(alpha=0.5, color='grey50',bin = 12, aes(y = after_stat(count / nrow(prova))))+labs(x='Patient mean switch percentage', y='Counts', title='Switch percentage distribution') + scale_fill_manual(name='',labels=c('Fast','Delayed'), values=c('#98d3c1','#eeb78c'))
#aes(y=after_stat(count/sum(count)))



prova<-prova[order(prova$sim_perc),,drop=F]
# divide at 55

qual_colors<-brewer.pal.info[brewer.pal.info$category=='qual',]
col_vector<-unlist(mapply(brewer.pal, qual_colors$maxcolors, rownames(qual_colors)))  
colors<-sample(col_vector, 2) 



surv_data <- survival
surv_data$delayed_switch_time <- as.integer(surv_data$clust)
surv_data[which(surv_data$delayed_switch_time == 2),]$delayed_switch_time <- 0

colnames(surv_data)[4] = "age_at_diagnosis"
param = "OS"
param.time = "OS.time"
fit = survfit(Surv(get(param.time), get(param)) ~ clust, data = surv_data)
plot(fit,xlab="Time (days)",ylab="Overall Survival",col=c('#eeb78c','#98d3c1'),lwd=2,cex=1.5)
stats <- coxph(Surv(get(param.time), get(param)) ~ delayed_switch_time + age_at_diagnosis + subtype, data = surv_data)
summary(stats)
ggforest(stats)

ggsurvplot(fit,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           legend.labs=c("delayed", "fast"), # change group labels
           legend.title="Average switch time",  # add legend title
           palette=c("#eeb78c", "#98d3c1"), # change colors of the groups
           title="Kaplan-Meier Curve for Breast Cancer Survival", # add title to plot
           risk.table.height=.2)
dev.off()



length(table(surv_data$age_at_diagnosis)) # 65

surv_data_switch <- prova
surv_data_switch$low_switch_commitment_percentage <- as.integer(surv_data_switch$clust2)
surv_data_switch[which(surv_data_switch$clust2 == 2),]$low_switch_commitment_percentage <- 1
surv_data_switch[which(surv_data_switch$clust2 == 1),]$low_switch_commitment_percentage <- 0

colnames(surv_data_switch)[4] = "age_at_diagnosis"
param = "OS"
param.time = "OS.time"
fit2 = survfit(Surv(get(param.time), get(param)) ~ clust2, data = surv_data_switch)
plot(fit2,xlab="Time (days)",ylab="Overall Survival",col=c("#98d3c1","#eeb78c"),lwd=2,cex=1.5)
stats2 <- coxph(Surv(get(param.time), get(param)) ~ low_switch_commitment_percentage + age_at_diagnosis + subtype, data = surv_data_switch)
summary(stats2)
ggforest(stats2)

ggsurvplot(fit2,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           legend.labs=c("high", "low"), # change group labels
           legend.title="Switch percentage",  # add legend title
           palette=c("#98d3c1", "#eeb78c"), # change colors of the groups
           title="Kaplan-Meier Curve for Breast Cancer Survival", # add title to plot
           risk.table.height=.2)
dev.off()
