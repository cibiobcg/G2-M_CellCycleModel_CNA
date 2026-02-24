library(ggplot2)
library(rstudioapi)
library(dplyr)
library(parallel)
library(data.table)
library(doParallel)
library(parallel)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


copy_numbers <- read.csv('../Files_generated/final_rates_copynumbers.csv')

copy_numbers$values <- round(copy_numbers$values,0)

system('mkdir Normal')
setwd('Normal')
normal <- data.frame()
#
function_normal <- function(val){
  system(paste0('../../Blenx_and_BetaWorkbench/sim64 ../../Blenx_and_BetaWorkbench/Cell-cycle_model.prog ../../Blenx_and_BetaWorkbench/Cell_cycle.types ../../Blenx_and_BetaWorkbench/Cell_cycle.fun -o=',as.factor(val),'_model_results'))

  upload <- fread(paste0(as.factor(val),'_model_results.E.out'), col.names =
                    c('time','cyclinB1_CDK1_p','cyclinB1_CDK1','CDC25C_p','CDC25C','WEE1_p','WEE1','FOXM1_plus_CTAK1','FCP1',
                      'C_TAK1','CHK1','Chp_i','Pho','null'), data.table = F,sep = '\t',fill = T)


  system(paste0('rm ',as.factor(val),'_model_results*'))
  if (any(as.numeric(upload$cyclinB1_CDK1) >= 150)){
    return(data.frame(Time = upload$time[which(upload$cyclinB1_CDK1 >= 150)][1], Specie = 'Normal'))
  }

}

normal <- rbind(normal,do.call(rbind, mclapply(X=1:100,FUN = function(X)function_normal(X),mc.cores = 10)))


setwd('..')

bef <- readLines('../Blenx_and_BetaWorkbench/Cell-cycle_model.prog')

Initial_val <- data.frame()

# #inhibitor_p part which is WEE1

species_fun <- function(iteration_x,value,specie,stat){
system(paste0('../../Blenx_and_BetaWorkbench/sim64 Model_Cell-cycle_',specie,'.prog ../../Blenx_and_BetaWorkbench/Cell_cycle.types ../../Blenx_and_BetaWorkbench/Cell_cycle.fun -o=',as.factor(iteration_x),'_Model_results_',specie,'_',as.factor(value)))

  upload <- fread(paste0(as.factor(iteration_x),'_Model_results_',specie,'_',as.factor(value),'.E.out'), col.names = c('time','cyclinB1_CDK1_p','cyclinB1_CDK1','CDC25C_p','CDC25C','WEE1_p','WEE1','FOXM1_plus_CTAK1','FCP1',
                                                                                                                       'C_TAK1','CHK1','Chp_i','Pho','null'), data.table = F,sep = '\t',fill = T)


  system(paste0('rm ',as.factor(iteration_x),'_Model_results_',specie,'_',as.factor(value),'*'))
  if (any(as.numeric(upload$cyclinB1_CDK1) >= 150)){
    print(data.frame(Time = upload$time[which(upload$cyclinB1_CDK1 >= 150)][1], Quant_spec = value, Specie = specie, Copy_number = stat))
    return(data.frame(Time = upload$time[which(upload$cyclinB1_CDK1 >= 150)][1], Quant_spec = value, Specie = specie, Copy_number = stat))
  }

}

# #
system('mkdir WEE1')
setwd('WEE1')

WEE1<- copy_numbers[which(copy_numbers$Gene == 'WEE1'),]

for (value in WEE1$values){
  
  status <- WEE1[which(WEE1$values == value),]$Copy_number
  word <- paste(as.character(value), 'WEE1', sep = ' ')
  print('Created the word')
  #print()

  bef1 <- gsub(pattern = "500 WEE1", replace = word, x = bef)
  print('Created the substitute')
  #
  writeLines(bef1, con = paste0('Model_Cell-cycle_WEE1.prog'))
  print('Writed new file')


  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'WEE1',status),mc.cores = 10)))

}

setwd('..')
system('mkdir FCP1')
setwd('FCP1')

CTPD1 <- copy_numbers[which(copy_numbers$Gene == "FCP1"),]
for (value in CTPD1$values) {
  
  status <- CTPD1[which(CTPD1$values == value),]$Copy_number
  word <- paste(as.character(value), 'FCP1', sep = ' ')
  bef1 <- gsub(pattern = "500 FCP1", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_FCP1.prog"))

  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'FCP1',status),mc.cores = 10)))
}

setwd('..')
system('mkdir C_TAK1')
setwd('C_TAK1')


MARK3 <- copy_numbers[which(copy_numbers$Gene == 'C-TAK1'),]
#
for(value in MARK3$values){
  
  status <- MARK3[which(MARK3$values == value),]$Copy_number
  word <- paste(as.character(value), 'C_TAK1', sep = ' ')
  bef1 <- gsub(pattern = "500 C_TAK1", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_C_TAK1.prog"))

  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'C_TAK1',status),mc.cores = 10)))

}

setwd('..')
system('mkdir FOXM1_plus_CTAK1')
setwd('FOXM1_plus_CTAK1')
#
FOXM1_plus_CTAK1 <- copy_numbers[which(copy_numbers$Gene == 'FOXM1' | copy_numbers$Gene == 'C-TAK1'),]
#
for(cp in unique(FOXM1_plus_CTAK1$Copy_number)){
  # print(cp)
  momentFOXM1 <- FOXM1_plus_CTAK1[which(FOXM1_plus_CTAK1$Gene == 'FOXM1' & FOXM1_plus_CTAK1$Copy_number == cp ),]
  momentCTAK1 <- FOXM1_plus_CTAK1[which(FOXM1_plus_CTAK1$Gene == 'C-TAK1' & FOXM1_plus_CTAK1$Copy_number == cp ),]
  if( cp == 'Shallow Deletion' | cp == "Diploid" ){
    if(momentCTAK1$values <= momentFOXM1$values){
      FOXM1_plus_CTAK1[which(FOXM1_plus_CTAK1$Gene == 'FOXM1' & FOXM1_plus_CTAK1$Copy_number == cp ),] <- NA
      # print('yep')
    }else{
      FOXM1_plus_CTAK1[which(FOXM1_plus_CTAK1$Gene == 'C-TAK1' & FOXM1_plus_CTAK1$Copy_number == cp ),] <- NA
      }
  }  else if(cp == 'Gain' | cp == 'Amplification'){
    if(momentCTAK1$values <= momentFOXM1$values){
      FOXM1_plus_CTAK1[which(FOXM1_plus_CTAK1$Gene == 'FOXM1' & FOXM1_plus_CTAK1$Copy_number == cp ),] <- NA
    }else{
      FOXM1_plus_CTAK1[which(FOXM1_plus_CTAK1$Gene == 'C-TAK1' & FOXM1_plus_CTAK1$Copy_number == cp ),] <- NA
    }
  }

}
#
final_FOXM1_plus_CTAK1<- FOXM1_plus_CTAK1 %>%  na.omit()
#
for (value in final_FOXM1_plus_CTAK1$values) {
  
  status <- final_FOXM1_plus_CTAK1[which(final_FOXM1_plus_CTAK1$values == value),]$Copy_number
  word <- paste(as.character(value), 'FOXM1_plus_CTAK1', sep = ' ')
  bef1 <- gsub(pattern = "500 FOXM1_plus_CTAK1", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_FOXM1_plus_CTAK1.prog"))

  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'FOXM1_plus_CTAK1',status),mc.cores = 10)))
  }


setwd('..')
system('mkdir CHK1')
setwd('CHK1')

CHK1 <- copy_numbers[which(copy_numbers$Gene == 'CHK1'),]

for (value in CHK1$values){

  status <- CHK1[which(CHK1$values == value),]$Copy_number
  word <- paste(as.character(value), 'CHK1', sep = ' ')
  bef1 <- gsub(pattern = "500 CHK1", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_CHK1.prog"))

  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'CHK1',status),mc.cores = 10)))

}

#parameters about kms 
setwd('..')
system('mkdir kms_cdk1_plus_CCNB1')
setwd('kms_cdk1_plus_CCNB1')


file_fun <- readLines('../../Blenx_and_BetaWorkbench/Cell_cycle.fun')

param_func <- function(iteration_x,val,parameter, stat){

  system(paste0('../../Blenx_and_BetaWorkbench/sim64 ../../Blenx_and_BetaWorkbench/Cell-cycle_model.prog ../../Blenx_and_BetaWorkbench/Cell_cycle.types Model_Cell-cycle.fun -o=',as.character(iteration_x),'_Model_results_',parameter,'_',as.character(val)))

  upload <- fread(paste0(as.character(iteration_x),'_Model_results_',parameter,'_',as.character(val),'.E.out'), col.names =c('time','cyclinB1_CDK1_p','cyclinB1_CDK1','CDC25C_p','CDC25C',
                      'WEE1_p','WEE1','FOXM1_plus_CTAK1','FCP1','C_TAK1','CHK1','Chp_i','Pho','null'), data.table = F,sep = '\t',fill = T)


  system(paste0('rm ',as.character(iteration_x),'_Model_results_',parameter,'_',as.character(val),'*'))
  if (any(as.numeric(upload$cyclinB1_CDK1) >= 150)){
    return(data.frame(Time = upload$time[which(upload$cyclinB1_CDK1 >= 150)][1], Quant_spec = val, Specie = parameter,Copy_number = stat))
  }


}
#
#
kms_cdk1_plus_CCNB1 <- copy_numbers[which(copy_numbers$Gene == 'CDK1'),]

for (rate in kms_cdk1_plus_CCNB1$linear_reg_values){

  status <- kms_cdk1_plus_CCNB1[which(kms_cdk1_plus_CCNB1$linear_reg_values == rate),]$Copy_number
  change <- (0.004/1.352279 ) * rate
  word <- paste0('let kms : const = ',as.factor(change),' ;')
  bef1 <- gsub(pattern = 'let kms : const = 0.004 ;', replace = word, x = file_fun)
  writeLines(bef1, con=paste0("Model_Cell-cycle.fun"))

  change <- round(change,5)
  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)param_func(X,change,'kms_cdk1_plus_CCNB1',status),mc.cores = 10)))
  

}

setwd('..')

#### Plot part 

library(ggplot2)
library(tidyverse)
library(paletteer)

`%nin%` <- Negate(`%in%`)

OskarSchlemmer <- paletteer::paletteer_d("lisa::OskarSchlemmer")

coloring<- c('Amplification' = '#3A488AFF','Gain' = '#8785B2FF','Diploid' = '#DABD61FF','Shallow Deletion' = '#D95F30FF','Deep Deletion' = '#BE3428FF')

Initial_val[which(Initial_val$Specie == 'kms_cdk1_plus_CCNB1'),]$Specie <- 'CCNB1-CDK1'
Initial_val_conserved <- Initial_val

`%nin%` <- Negate(`%in%`)

ggplot(data=Initial_val %>% group_by(Specie, Quant_spec,Copy_number) %>% summarise(Mean = mean(Time)), aes(Mean, Specie)) +
  geom_boxplot() + geom_point(size=3, aes(Mean,Specie,colour = Copy_number)) +scale_color_manual(values = coloring)+
  geom_vline(xintercept = mean(normal$Time), colour = '#b91000', lwd =0.4) + 
  labs(title = 'Copy number model switching times',x='Time' ,y='Species',caption = paste0('Normal status:',as.character(round(mean(normal$Time),1))), colour = 'Status') 

success_perc <- Initial_val %>% group_by(Specie, Quant_spec,Copy_number) %>% dplyr::count(Quant_spec) %>% mutate(perc = n)

level_order <- c('Deep Deletion', 'Shallow Deletion', 'Diploid','Gain','Amplification') 

ggplot(success_perc, aes(factor(Copy_number,levels = level_order), Specie, fill= perc)) + 
  geom_tile(color = 'black') +
  geom_text(aes(label = perc), color = "white", size = 2.5) +
  scale_fill_gradient(low = '#605b8f',high = '#de425b') +
  guides(fill = guide_colorbar(title = 'Percentage')) +
  labs(title = 'Copy number model percentage of switching',caption = paste0('Normal status:',as.character(nrow(normal)),'%'),x = 'Copy Number',y = "Species")
