library(ggplot2)
library(rstudioapi)
library(dplyr)
library(parallel)
library(data.table)
library(doParallel)
library(parallel)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

system('mkdir normal')
setwd('normal/')
`%nin%` <- Negate(`%in%`)


normal <- data.frame()


function_normal <- function(val,TF){

  system(paste0('../../Blenx_and_BetaWorkbench/sim64 ../../Blenx_and_BetaWorkbench/Cell-cycle_model.prog ../../Blenx_and_BetaWorkbench/Cell_cycle.types ../../Blenx_and_BetaWorkbench/Cell_cycle.fun -o=',as.factor(val),'_model_results'))
  
  upload <- fread(paste0(as.factor(val),'_model_results.E.out'), col.names =
                    c('time','cyclinB1_CDK1_p','cyclinB1_CDK1','CDC25C_p','CDC25C','WEE1_p','WEE1','FOXM1_plus_CTAK1','FCP1',
                      'C_TAK1','CHK1','Chp_i','Pho','null'), data.table = F,sep = '\t',fill = T)
  

  system(paste0('rm ',as.factor(val),'_model_results_',TF,'*'))
  if (any(as.numeric(upload$TR_p) >= 150)){
    return(data.frame(Time = upload$time[which(upload$TR_p >= 150)][1], Specie = 'Normal', TF_a = TF))
  }
  
}

normal <- rbind(normal,do.call(rbind, mclapply(X=1:100,FUN = function(X)function_normal(X,'500'),mc.cores = 10)))


ggplot(normal, aes(Time)) + geom_histogram(binwidth = 50,color = "#000000", fill = "#0099F8") + labs(title = 'Switch Time Distribution')

setwd('..')
# 
bef <- readLines('../Blenx_and_BetaWorkbench/Cell-cycle_model.prog')

Initial_val <- data.frame()


species_fun <- function(iteration_x,value,specie,TF){
  system(paste0('../../Blenx_and_BetaWorkbench/sim64 Model_Cell-cycle_',TF,'_',specie,'.prog../../Blenx_and_BetaWorkbench/Cell_cycle.types ../../Blenx_and_BetaWorkbench/Cell_cycle.fun -o=',as.factor(iteration_x),'_Model_results_',TF,'_',specie,'_',as.factor(value)))
  
  upload <- fread(paste0(as.factor(iteration_x),'_Model_results_',specie,'_',as.factor(value),'.E.out'), col.names = c('time','cyclinB1_CDK1_p','cyclinB1_CDK1','CDC25C_p','CDC25C','WEE1_p','WEE1','FOXM1_plus_CTAK1','FCP1',
                                                                                                                       'C_TAK1','CHK1','Chp_i','Pho','null'), data.table = F,sep = '\t',fill = T)
  
  
  system(paste0('rm ',as.factor(iteration_x),'_Model_results_',TF,'_',specie,'_',as.factor(value),'*'))
  if (any(as.numeric(upload$cyclinB1_CDK1) >= 150)){
    print(data.frame(Time = upload$time[which(upload$cyclinB1_CDK1 >= 150)][1], Quant_spec = value, Specie = specie, TF_a = TF))
    return(data.frame(Time = upload$time[which(upload$cyclinB1_CDK1 >= 150)][1], Quant_spec = value, Specie = specie, TF_a = TF))
  }
  
}


system('mkdir WEE1')
setwd('WEE1')

for (value in seq(0,1000,50)){
    word <- paste(as.factor(value), 'WEE1', sep = ' ')

    bef1 <- gsub(pattern = "500 WEE1", replace = word, x = bef)
    writeLines(bef1, con = paste0("Model_Cell-cycle_500_WEE1.prog"))

    Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'WEE1','500'),mc.cores = 10)))
    
}

setwd('..')
system('mkdir FCP1')
setwd('FCP1')
# 
for (value in seq(0,1000,50)) {
  word <- paste(as.factor(value), 'FCP1', sep = ' ')

  bef1 <- gsub(pattern = "500 FCP1", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_500_FCP1.prog"))

  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'FCP1','500'),mc.cores = 10)))
}

setwd('..')
system('mkdir C_TAK1')
setwd('C_TAK1')
# 
for(value in seq(0,1000,50)){
  word <- paste(as.factor(value), 'C_TAK1', sep = ' ')
  bef1 <- gsub(pattern = "500 C_TAK1", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_500_C_TAK1.prog"))

  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'C_TAK1','500'),mc.cores = 10)))
  
}

setwd('..')
system('mkdir FOXM1_plus_CTAK1')
setwd('FOXM1_plus_CTAK1')

for (value in seq(0,1000,50)) {
  word <- paste(as.factor(value), 'FOXM1_plus_CTAK1', sep = ' ')

  bef1 <- gsub(pattern = "500 FOXM1_plus_CTAK1", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_500_FOXM1_plus_CTAK1.prog"))

  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'FOXM1_plus_CTAK1','500'),mc.cores = 10)))
}


setwd('..')
system('mkdir Chp_i')
setwd('Chp_i')

for (value in seq(0,1000,50)){
  word <- paste(as.factor(value), 'Chp_i', sep = ' ')
  
  bef1 <- gsub(pattern = "500 Chp_i", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_500_Chp_i.prog"))
  
  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'Chp_i','500'),mc.cores = 10)))
  
}

setwd('..')
system('mkdir Pho')
setwd('Pho')

for (value in seq(0,1000,50)){
  word <- paste(as.factor(value), 'Pho', sep = ' ')
  
  bef1 <- gsub(pattern = "500 Pho", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_500_Pho.prog"))
  
  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'Pho','500'),mc.cores = 10)))
  
}


setwd('..')
system('mkdir CHK1')
setwd('CHK1')

for (value in seq(0,1000,50)){
  word <- paste(as.factor(value), 'CHK1', sep = ' ')
  
  bef1 <- gsub(pattern = "500 CHK1", replace = word, x = bef)
  writeLines(bef1, con = paste0("Model_Cell-cycle_500_CHK1.prog"))
  
  Initial_val <- rbind(Initial_val, do.call(rbind,mclapply(X=1:100,function(X)species_fun(X,value,'CHK1','500'),mc.cores = 10)))
  
}

#parameters
setwd('..')
system('mkdir parameters')
setwd('parameters')
###
par_500 <- data.frame()

file_fun <- readLines('../../Blenx_and_BetaWorkbench/Cell_cycle.fun')

param_func <- function(iteration_x,val,parameter,TF){
  
  system(paste0('../../Blenx_and_BetaWorkbench/sim64 ../../Blenx_and_BetaWorkbench/Cell-cycle_model.prog ../../../Cell_cycle.types Model_Cell-cycle_',TF,'.fun -o=',as.factor(iteration_x),'_model_results_',TF,'_',parameter,'_',as.factor(val)))
  
 
  upload <- fread(paste0(as.factor(iteration_x),'_model_results_',TF,'_',parameter,'_',as.factor(value),'.E.out'), col.names = 
                    c('time','cyclinB1_CDK1_p','cyclinB1_CDK1','CDC25C_p','CDC25C','WEE1_p','WEE1','FOXM1_plus_CTAK1','FCP1',
                        'C_TAK1','CHK1','Chp_i','Pho','null'), data.table = F,sep = '\t',fill = T)
  
  system(paste0('rm ',as.factor(iteration_x),'_model_results_',TF,'_',parameter,'_',as.factor(val),'*'))
  if (any(as.numeric(upload$TR_p) >= 150)){
    return(data.frame(Time = upload$time[which(upload$TR_p >= 150)][1], Quant_par = val, Parameter = parameter, TF_a = TF))
  }
  

}


for (par in file_fun[1:40]){
  value = strsplit(par,' ')[[1]][6]
  parameter = strsplit(par,' ')[[1]][2]
  if (value != '1/500' ){
    value = as.numeric(value)
    for (j in c(seq(0.1,1,0.1),seq(2,10,1))){
      change = j * value
      word <- paste0('let ',parameter,' : const = ',as.factor(change),' ;')
      bef1 <- gsub(pattern = par, replace = word, x = file_fun)
      writeLines(bef1, con=paste0("Model_Cell-cycle_500.fun"))
      
      par_500 <- rbind(par_500, do.call(rbind,mclapply(X=1:100,function(X)param_func(X,change,parameter,'500'),mc.cores = 10)))
      
      }
  }
  
}

# plotting species #####

#times
ggplot(data=Initial_val %>% group_by(Specie, Quant_spec) %>% summarise(Mean = mean(Time)), aes(Mean, Specie, fill = Specie)) +
  geom_boxplot() + geom_vline(xintercept = mean(normal$Time), colour = 'red', lwd =0.2) +
  geom_point(position=position_dodge(),aes(group = Mean, colour = Quant_spec)) +
  scale_colour_gradient(low = '#605b8f',high = '#de425b') + labs(title = 'Sensitivity Analysis Species',x='Time' ,y='Species')


moment <- Initial_val %>% group_by(Specie, Quant_spec) %>% dplyr::count(Quant_spec) %>% mutate(perc = n)
missing_data <- data.frame(Specie = NA,Quant_spec = NA,n = NA, perc = NA)
for (spec in unique(moment$Specie)){
  selected_spec = moment[which(moment$Specie == spec),]$Quant_spec
  for (val in seq(0,1000,50)){
    if (val %nin% selected_spec){
      mom <- data.frame(Specie = spec, Quant_spec = val, n = 0 , perc = 0 )
      missing_data <- rbind(missing_data,mom)
    }
  }
}

moment <- rbind(moment,missing_data %>% na.omit())

ggplot(moment, aes(Quant_spec, Specie, fill= perc)) + 
  geom_tile(color = 'black') +
  geom_text(aes(label = perc), color = "white", size = 2.5) +
  scale_fill_gradient(low = '#605b8f',high = '#de425b') +
  guides(fill = guide_colorbar(title = 'Percentage')) +
  labs(title = 'Sensitivity Analysis Species',caption = paste0('Normal status:',as.character(nrow(normal)),'%'),x = 'Quantity specie',y = "Species")



# plotting parameters  #####

# time
parameters_of_interes <- c('kms','kma','jma','kmi','jmi','kmd','kmd_1','jwa','kwi','perc','jwi','kca','jca','kci','jci','kcd','kcs','kcd_1','kcp_7','jcp_7','kcp_6','jcp_6','jcp_1',
                           'kcp_1','kcp_2','jcp_2','jcp_5','kcp_5','jcp_3','kcp_3','jcp_4','kcp_4','kmi_1','kcp_8','jcp_8')

par_500 <- par_500[which(par_500$Parameter %in% parameters_of_interes),]

par_500 <- par_500 %>% mutate(Slots = case_when(Parameter == 'kms' ~ Quant_par / 0.004,
                                                Parameter == 'kma' ~ Quant_par / 0.5,
                                                Parameter == 'jma' ~ Quant_par / 1,
                                                Parameter == 'kmi' ~ Quant_par  / 0.5,
                                                Parameter == 'jmi' ~ Quant_par  / 1,
                                                Parameter == 'kmd' ~ Quant_par  / 0.002,
                                                Parameter == 'kmd_1' ~ Quant_par / 0.002,
                                                Parameter == 'jwa' ~ Quant_par / 0.1,
                                                Parameter == 'kwi' ~ Quant_par / 2,
                                                Parameter == 'perc' ~ Quant_par / 0.01,
                                                Parameter == 'jwi' ~ Quant_par / 0.1,
                                                Parameter == 'kca' ~ Quant_par / 2,
                                                Parameter == 'jca' ~ Quant_par / 0.1,
                                                Parameter == 'kci' ~ Quant_par / 0.2,
                                                Parameter == 'jci' ~ Quant_par / 0.1,
                                                Parameter == 'kcd' ~ Quant_par / 0.01,
                                                Parameter == 'kcs' ~ Quant_par / 0.01,
                                                Parameter == 'kcd_1' ~ Quant_par / 0.01,
                                                Parameter == 'kcp_7' ~ Quant_par / 0.2,
                                                Parameter == 'jcp_7' ~ Quant_par / 0.1,
                                                Parameter == 'kcp_6' ~ Quant_par / 0.2,
                                                Parameter == 'jcp_6' ~ Quant_par / 0.1,
                                                Parameter == 'kcp_1' ~ Quant_par / 0.1,
                                                Parameter == 'jcp_1' ~ Quant_par / 2,
                                                Parameter == 'kcp_2' ~ Quant_par / 0.2,
                                                Parameter == 'jcp_2' ~ Quant_par / 0.1,
                                                Parameter == 'kcp_5' ~ Quant_par / 0.2,
                                                Parameter == 'jcp_5' ~ Quant_par / 0.1,
                                                Parameter == 'kcp_4' ~ Quant_par / 2,
                                                Parameter == 'jcp_4' ~ Quant_par / 0.1,
                                                Parameter == 'kcp_3' ~ Quant_par / 2,
                                                Parameter == 'jcp_3' ~ Quant_par / 0.1,
                                                Parameter == 'kmi_1' ~ Quant_par / 2,
                                                Parameter == 'kcp_8' ~ Quant_par / 0.1,
                                                Parameter == 'jcp_8' ~ Quant_par / 0.1,
))

ggplot(data=par_500 %>% group_by(Parameter, Quant_par,Slots) %>% summarise(Mean = mean(Time)), aes(Mean, Parameter, fill = Parameter)) +
  geom_boxplot() + geom_point(position=position_dodge(),aes(group = Mean, colour = Slots))  +
  geom_vline(xintercept = mean(normal500$Time), colour = '#b91000', lwd =0.4) + 
  scale_colour_gradient(low = '#605b8f',high = '#de425b') + labs(title = 'Sensitivity Analysis Parameters',x='Time' ,y='Parameters',caption = paste0('Normal status:',as.character(round(mean(normal$Time),2)))) 


# percentage 

moment1 <- par_500 %>% group_by(Parameter, Quant_par,Slots) %>% dplyr::count(Quant_par) %>% mutate(perc = n)

missing_data <- data.frame(Parameter = NA,Quant_par = NA,n = NA, perc = NA, Slots =NA)
for (par in unique(moment1$Parameter)){
  selected_par = as.character(moment1[which(moment1$Parameter == par),]$Slots)
  for (val in c(seq(0.1,1,0.1),seq(2,10,1))){
    if (as.character(val) %nin% selected_par){
      mom <- data.frame(Parameter = par, Quant_par = 0, n = 0 , perc = 0, Slots = val)
      missing_data <- rbind(missing_data,mom)
      
    }
  }
}

moment1 <- rbind(moment1,missing_data %>% na.omit())
moment1$Slots <- as.character(moment1$Slots)

ggplot(moment1, aes(Slots, Parameter, fill= perc)) + 
  geom_tile(color = 'black') +
  geom_text(aes(label = perc), color = "white", size = 3) +
  scale_fill_gradient(low = '#605b8f',high = '#de425b') +
  scale_x_discrete(limits = c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9', '1','2','3','4','5','6','7','8','9','10')) +
  guides(fill = guide_colorbar(title = 'Percentage')) +
  labs(title = 'Sensitivity Analysis Parameters',caption = paste0('Normal status:',as.character(nrow(normal)),'%'),x = 'Slots',y = "Parameters")

