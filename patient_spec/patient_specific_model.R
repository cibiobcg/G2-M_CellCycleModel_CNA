library(ggplot2)
library(rstudioapi)
library(dplyr)
library(parallel)
library(data.table)
library(doParallel)
library(parallel)
library(xlsx)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

patient_spec <- read.csv('../Files_generated/table_patient_specific.csv')

copy_number <- read.csv('../Files_generated/final_rates_copynumbers.csv')

bef <- readLines('../Blenx_and_BetaWorkbench/Cell-cycle_model.prog')
file_fun <- readLines('../Blenx_and_BetaWorkbench/Cell_cycle.fun')

system('mkdir patient_folder')
setwd('patient_folder')

Initial_val <- data.frame()


function_patient <- function(val, patient,notes){
  system(paste0('../../Blenx_and_BetaWorkbench/sim64 Cell-cycle_model_changed.prog ../../Blenx_and_BetaWorkbench/Cell_cycle.types Cell_cycle_changed.fun -o=',as.factor(val),'_model_results_',patient))
  
  upload <- fread(paste0(as.factor(val),'_model_results_',patient,'.E.out'), col.names =
                    c('time','cyclinB1_CDK1_p','cyclinB1_CDK1','CDC25C_p','CDC25C','WEE1_p','WEE1','FOXM1_plus_CTAK1','FCP1',
                      'C_TAK1','CHK1','Chp_i','Pho','null'), data.table = F,sep = '\t',fill = T)
  
  
  system(paste0('rm ',as.factor(val),'_model_results_',patient,'*'))
  if (any(as.numeric(upload$cyclinB1_CDK1) >= 150)){
    return(data.frame(Time = upload$time[which(upload$cyclinB1_CDK1 >= 150)][1], Patient = patient, Notes_situation_cp = notes))
  }
  
}

counter = 0 
for (patient in unique(patient_spec$sample)){
  counter = counter + 1
  
  print(paste0('Patient number ',counter,' out of 983'))
  momentary <- patient_spec[which(patient_spec$sample == patient),]

  status_wee1 <- momentary[which(momentary$hgnc_symbol == 'WEE1'),]$copy_number_alteration
  status_foxm1 <- momentary[which(momentary$hgnc_symbol == 'FOXM1'),]$copy_number_alteration
  status_ctak1 <- momentary[which(momentary$hgnc_symbol == 'C-TAK1'),]$copy_number_alteration
  status_fcp1 <- momentary[which(momentary$hgnc_symbol == 'FCP1'),]$copy_number_alteration
  status_chk1 <- momentary[which(momentary$hgnc_symbol == 'CHK1'),]$copy_number_alteration
  status_cdk1 <- momentary[which(momentary$hgnc_symbol == 'CDK1'),]$copy_number_alteration
  
  # i take the one that present less transcript 
  if(copy_number[which(copy_number$Gene == 'FOXM1' & copy_number$Copy_number == status_foxm1),]$value <= copy_number[which(copy_number$Gene == 'C-TAK1' & copy_number$Copy_number == status_ctak1),]$value){
    val_FOXM1_plus_CTAK1 <- copy_number[which(copy_number$Gene == 'FOXM1' & copy_number$Copy_number == status_foxm1),]$value
  }else{
    
    val_FOXM1_plus_CTAK1 <- copy_number[which(copy_number$Gene == 'C-TAK1' & copy_number$Copy_number == status_ctak1),]$value
  }
  
  new_last_line <- paste0('run ', as.character(round(copy_number[which(copy_number$Gene == 'WEE1' & copy_number$Copy_number == status_wee1),]$value,0)),' WEE1 || ',as.character(round(val_FOXM1_plus_CTAK1,0)),' FOXM1_plus_CTAK1 || ',
                          as.character(round(copy_number[which(copy_number$Gene == 'FCP1' & copy_number$Copy_number == status_fcp1),]$value,0)),' FCP1 || ',as.character(round(copy_number[which(copy_number$Gene == 'C-TAK1' & copy_number$Copy_number == status_ctak1),]$value,0)),
                          ' C_TAK1 || ',as.character(round(copy_number[which(copy_number$Gene == 'CHK1' & copy_number$Copy_number == status_chk1),]$value,0)),' CHK1 || 500 Chp_i || 500 Pho  ')
  
  bef1 <- gsub(pattern = "run 500 WEE1 \\|\\| 500 FOXM1_plus_CTAK1 \\|\\| 500 FCP1 \\|\\| 500 C_TAK1 \\|\\| 500 CHK1 \\|\\| 500 Chp_i \\|\\| 500 Pho  ", replacement = new_last_line, x = bef)
  
  writeLines(bef1, con = paste0("Cell-cycle_model_changed.prog"))
  
  rate_kms <- copy_number[which(copy_number$Gene == 'CDK1' & copy_number$Copy_number == status_cdk1),]$linear_reg_values
  divider_kms <- copy_number[which(copy_number$Gene == 'CDK1' & copy_number$Copy_number == 'Diploid'),]$linear_reg_values
  
  change <- (0.004/divider_kms) * rate_kms
  word_kms <- paste0('let kms : const = ',as.factor(change),' ;')
  bef_kms <- gsub(pattern = 'let kms : const = 0.004 ;', replace = word_kms, x = file_fun)
  writeLines(bef_kms, con=paste0("Cell_cycle_changed.fun"))

  note_line <- paste0('WEE1 ',status_wee1, ' | FOXM1 ', status_foxm1, ' | C-TAK1 ',status_ctak1,' | FCP1 ',status_fcp1,' | CHK1 ',status_chk1,' | CDK1 ',status_cdk1)
  # print()
  Initial_val <- rbind(Initial_val,do.call(rbind, mclapply(X=1:100,FUN = function(X)function_patient(X,patient,note_line),mc.cores = 10)))
  
}

## 

library(ggplot2)
library(tidyverse)
library(paletteer)

`%nin%` <- Negate(`%in%`)

metadata_patient <- read_tsv('../../Files_generated/brca_tcga_pan_can_atlas_2018_clinical_data.tsv')

patients_final <- Initial_val  

patient_med_time <- patients_final %>% group_by(Patient, Notes_situation_cp) %>% summarise(Mean = mean(Time))

patient_perc <- patients_final %>% group_by(Patient, Notes_situation_cp) %>% dplyr::count(Patient) %>% mutate(Percentage_switch = n)

patient_med_time$Percentage_switch <- patient_perc$Percentage_switch

patient_original <- patient_spec

not_presents <- patient_original[which(patient_original$sample %nin% patient_med_time$Patient),]

keep_missing <- data.frame()

for (pat in unique(not_presents$sample)){
  
  momentary <- not_presents[which(not_presents$sample == pat),]
  
  notes_info<-''
  
  for(genes in momentary$hgnc_symbol){
    status <- momentary[which(momentary$hgnc_symbol == genes),]$copy_number_alteration
    
    notes_info <- paste0(notes_info,' ',genes,' ',status, '|')
  }
  
  moment2 <- data.frame(Patient = pat,Notes_situation_cp = notes_info, Mean = NA , Percentage_switch = 0 )
  
  keep_missing <- rbind(keep_missing,moment2)
}


final_patients <- rbind(patient_med_time,keep_missing)

metadata_filtered <- metadata_patient[which(metadata_patient$`Sample ID` %in% final_patients$Patient),]

write.csv(final_patients,'../../Files_generated/final_patient.csv')

