library(org.Hs.eg.db)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(recount3)
library(recount)
library(biomaRt)
library(sva)
library(rstudioapi)
library(ggplot2)
library(preprocessCore)
library(ggplot2)
library(dplyr)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
gene_names <- list('FOXM1','CDC25C','CCNB1','CHEK1','WEE1','CDK1','MARK3','CTDP1') 
`%nin%` <- Negate(`%in%`)



brca1 <- recount3::create_rse_manual(
  project = "BRCA",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

assays(brca1)$counts <- compute_read_counts(brca1)

##### retrive tcga barcodes of the samples ####

tcga <- as.data.frame(brca1@colData)

tcga_normal <- tcga[which(tcga$tcga.cgc_sample_sample_type == 'Solid Tissue Normal'),]


tcga <- tcga[which(tcga$tcga.cgc_sample_sample_type != 'Solid Tissue Normal' & tcga$tcga.cgc_sample_sample_type != 'Metastatic' ),] # removed the normal samples and remove also the Metastatic 


#there are samples that present the same tcga barcod, retrive the rows 

dup <- tcga[duplicated(tcga$tcga.gdc_cases.submitter_id),]# there are duplicates 

# for normal tissues
new_barcodes_healthy <- list()
n = 0
for (bar in tcga_normal$tcga.tcga_barcode){
  if(grepl('11[[:alpha:]]',bar)){
    new_barcodes_healthy <- append(new_barcodes_healthy,bar)
  } }


newer_barcodes_healthy <- list()
for (barcode in new_barcodes_healthy){
  phases <- unlist(strsplit(barcode,'-',fixed = T))
  change_bar <- substr(phases[4],1,2)
  bar_new <- paste0(phases[1],'-',phases[2],'-',phases[3],'-',change_bar)
  newer_barcodes_healthy <- append(newer_barcodes_healthy,bar_new)
  
}

# tumor 
new_barcodes <- list()
n = 0
for (bar in tcga$tcga.tcga_barcode){
  if(grepl('01[[:alpha:]]',bar)){
    new_barcodes <- append(new_barcodes,bar)
  } }


#cambiare e tenere solo quelli con A finale
newer_barcodes <- list()
for (barcode in new_barcodes){
  phases <- unlist(strsplit(barcode,'-',fixed = T))
  change_bar <- substr(phases[4],1,3)
  bar_new <- paste0(phases[1],'-',phases[2],'-',phases[3],'-',change_bar)
  newer_barcodes <- append(newer_barcodes,bar_new)
  
}
###TPM ####
Brca_tpm <- as.data.frame(recount::getTPM(brca1,length_var = 'bp_length', mapped_var = NULL))

Brca_tpm2 <- str_replace(rownames(Brca_tpm),
                         pattern = ".[0-9]+$",
                         replacement = "")
rownames(Brca_tpm) <- Brca_tpm2

Brca_tpm_norm <- Brca_tpm[which(colnames(Brca_tpm) %in% tcga_normal$external_id)]

Brca_tpm <- Brca_tpm[which(colnames(Brca_tpm) %in% tcga$external_id)]

Brca_tpm <- Brca_tpm %>% rownames_to_column('ensembl_id')
Brca_tpm_norm <- Brca_tpm_norm %>%  rownames_to_column('ensembl_id')

colnames(Brca_tpm)[2:1128] <- newer_barcodes
colnames(Brca_tpm_norm)[2:113] <-newer_barcodes_healthy


### for tumors there are replicates so 
Brca_tpm <- Brca_tpm[, !duplicated(colnames(Brca_tpm))]

Brca_tpm_momentary<- Brca_tpm %>% dplyr::select(ends_with('A'))


newest_barcodes <- list()
for (barcode in colnames(Brca_tpm_momentary)){
  phases <- unlist(strsplit(barcode,'-',fixed = T))
  change_bar <- substr(phases[4],1,2)
  bar_new <- paste0(phases[1],'-',phases[2],'-',phases[3],'-',change_bar)
  newest_barcodes <- append(newest_barcodes,bar_new)
  
}

colnames(Brca_tpm_momentary) <-newest_barcodes

Brca_tpm_momentary$ensembl_id <- Brca_tpm$ensembl_id # 1080 samples

trial1 <- sapply(Brca_tpm_momentary[1:1080], as.matrix)
trail2 <- tidyr::gather(as.data.frame(trial1[,1:80]),key = 'sample',value = 'read_number')

ggplot() +
  geom_boxplot(alpha = 0.7, mapping = aes(sample,log10(read_number+1)),data = trail2)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))
  

### healthy tissues
trial1_norm <- sapply(Brca_tpm_norm[2:113], as.matrix)
trail2_norm <- tidyr::gather(as.data.frame(trial1_norm[,1:80]),key = 'sample',value = 'read_number')

ggplot() +
  geom_boxplot(alpha = 0.7, mapping = aes(sample,log10(read_number+1)),data = trail2_norm)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))

#### quantile normalization ####
Brca_norm <- as.data.frame(normalize.quantiles(sapply(Brca_tpm_momentary[1:1080], as.matrix), copy = TRUE))
colnames(Brca_norm) <- colnames(Brca_tpm_momentary[1:1080])
Brca_norm_plot <- tidyr::gather(Brca_norm[,1:80],key = 'sample',value = 'read_number')
 
ggplot() +
  geom_boxplot(alpha = 0.7, mapping = aes(sample,log10(read_number+1)),data = Brca_norm_plot)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))

### healthy tissues 
Brca_norm_healthy <- as.data.frame(normalize.quantiles(sapply(Brca_tpm_norm[2:113], as.matrix), copy = TRUE))
colnames(Brca_norm_healthy) <- colnames(Brca_tpm_norm[2:113])
Brca_norm_plot_healthy <- tidyr::gather(Brca_norm_healthy[,1:80],key = 'sample',value = 'read_number')
# 
ggplot() +
  geom_boxplot(alpha = 0.7, mapping = aes(sample,log10(read_number+1)),data = Brca_norm_plot_healthy)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))#+


Brca_norm <- Brca_norm %>% add_column(ensembl_gene_id = Brca_tpm_momentary$ensembl_id, .before = 'TCGA-E9-A249-01')
Brca_norm_healthy <- Brca_norm_healthy  %>% add_column(ensembl_gene_id = Brca_tpm_norm$ensembl_id, .before = 'TCGA-BH-A0B8-11')


ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'hsapiens_gene_ensembl',
                      version = 113)
annotation_brca <- getBM(attributes = c( 'hgnc_symbol','ensembl_gene_id'),
                         filters = 'hgnc_symbol',
                         values = unlist(gene_names), 
                         mart = ensembl)

Brca_final <- merge(Brca_norm,annotation_brca,by = 'ensembl_gene_id')
row_sums <- rowSums(Brca_final[,2:1080])
row_sums # second CTDP1 present less information so will be eliminated 
Brca_final <- Brca_final[!duplicated(Brca_final$hgnc_symbol),]

Brca_norm_plot <- tidyr::gather(Brca_final,key = 'sample',value = 'value',-c('ensembl_gene_id','hgnc_symbol'))
Brca_norm_plot$value <- as.numeric(Brca_norm_plot$value)
Brca_norm_plot$value <- round(Brca_norm_plot$value,2)

ggplot(Brca_norm_plot, aes(x=hgnc_symbol, y =log10(value +1) , fill =hgnc_symbol )) + geom_boxplot() + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Gene Name') + ylab('log10(TPM + 1)') + ggtitle('BRCA TCGA 1080 samples TPM then quantile normalization')

####### healthy 

Brca_final_healthy <- merge(Brca_norm_healthy,annotation_brca,by = 'ensembl_gene_id')
row_sums_healthy <- rowSums(Brca_final_healthy[,2:113])
Brca_final_healthy <- Brca_final_healthy[!duplicated(Brca_final_healthy$hgnc_symbol),]

Brca_norm_plot_healthy_2 <- tidyr::gather(Brca_final_healthy,key = 'sample',value = 'value',-c('ensembl_gene_id','hgnc_symbol'))
Brca_norm_plot_healthy_2$value <- round(Brca_norm_plot_healthy_2$value,2)

ggplot(Brca_norm_plot_healthy_2, aes(x=hgnc_symbol, y =log10(value +1) , fill =hgnc_symbol )) + geom_boxplot() + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Gene Name') + ylab('log10(TPM + 1)') + ggtitle('BRCA TCGA 112 healthy samples TPM then quantile normalization')

##########################

##########################
setwd('Cbioportal')

FOXM1 <- read.delim('FOXM1_cbiportal_dataplot.txt',header = T) 

FOXM1 <- FOXM1[!duplicated(FOXM1$Sample.Id), ]

Brca_norm_plot$copy_number_alteration <- NA

for (sample in FOXM1$Sample.Id){
  if(sample %in% Brca_norm_plot$sample){
    Brca_norm_plot[which(Brca_norm_plot$sample == sample & Brca_norm_plot$hgnc_symbol == 'FOXM1'),]$copy_number_alteration <- FOXM1[which(FOXM1$Sample.Id == sample),]$FOXM1..Putative.copy.number.alterations.from.GISTIC
     
  }
}

CCNB1 <- read.delim('CCNB1_cbiportal_dataplot.txt',header = T) 
CCNB1 <- CCNB1[!duplicated(CCNB1$Sample.Id), ]

for (sample in CCNB1$Sample.Id){
  if(sample %in% Brca_norm_plot$sample){
    Brca_norm_plot[which(Brca_norm_plot$sample == sample & Brca_norm_plot$hgnc_symbol == 'CCNB1'),]$copy_number_alteration <- CCNB1[which(CCNB1$Sample.Id == sample),]$CCNB1..Putative.copy.number.alterations.from.GISTIC
    
  }
}

CDC25C <- read.delim('CDC25C_cbiportal_dataplot.txt',header = T) 
CDC25C <- CDC25C[!duplicated(CDC25C$Sample.Id), ]

for (sample in CDC25C$Sample.Id){
  if(sample %in% Brca_norm_plot$sample){
    Brca_norm_plot[which(Brca_norm_plot$sample == sample & Brca_norm_plot$hgnc_symbol == 'CDC25C'),]$copy_number_alteration <- CDC25C[which(CDC25C$Sample.Id == sample),]$CDC25C..Putative.copy.number.alterations.from.GISTIC
    
  }
}

CDK1 <- read.delim('CDK1_cbiportal_dataplot.txt',header = T) 
CDK1 <- CDK1[!duplicated(CDK1$Sample.Id), ]

for (sample in CDK1$Sample.Id){
  if(sample %in% Brca_norm_plot$sample){
    Brca_norm_plot[which(Brca_norm_plot$sample == sample & Brca_norm_plot$hgnc_symbol == 'CDK1'),]$copy_number_alteration <- CDK1[which(CDK1$Sample.Id == sample),]$CDK1..Putative.copy.number.alterations.from.GISTIC
    
  }
}

CHEK1 <- read.delim('CHEK1_cbiportal_dataplot.txt',header = T) 
CHEK1 <- CHEK1[!duplicated(CHEK1$Sample.Id), ]

for (sample in CHEK1$Sample.Id){
  if(sample %in% Brca_norm_plot$sample){
    Brca_norm_plot[which(Brca_norm_plot$sample == sample & Brca_norm_plot$hgnc_symbol == 'CHEK1'),]$copy_number_alteration <- CHEK1[which(CHEK1$Sample.Id == sample),]$CHEK1..Putative.copy.number.alterations.from.GISTIC
    
  }
}

CTDP1 <- read.delim('CTDP1_cbiportal_dataplot.txt',header = T) 
CTDP1 <- CTDP1[!duplicated(CTDP1$Sample.Id), ]

for (sample in CTDP1$Sample.Id){
  if(sample %in% Brca_norm_plot$sample){
    Brca_norm_plot[which(Brca_norm_plot$sample == sample & Brca_norm_plot$hgnc_symbol == 'CTDP1'),]$copy_number_alteration <- CTDP1[which(CTDP1$Sample.Id == sample),]$CTDP1..Putative.copy.number.alterations.from.GISTIC
    
  }
}

MARK3 <- read.delim('MARK3_cbiportal_dataplot.txt',header = T) 
MARK3 <- MARK3[!duplicated(MARK3$Sample.Id), ]

for (sample in MARK3$Sample.Id){
  if(sample %in% Brca_norm_plot$sample){
    Brca_norm_plot[which(Brca_norm_plot$sample == sample & Brca_norm_plot$hgnc_symbol == 'MARK3'),]$copy_number_alteration <- MARK3[which(MARK3$Sample.Id == sample),]$MARK3..Putative.copy.number.alterations.from.GISTIC
    
  }
}

WEE1 <- read.delim('WEE1_cbiportal_dataplot.txt',header = T) 
WEE1 <- WEE1[!duplicated(WEE1$Sample.Id), ]

for (sample in WEE1$Sample.Id){
  if(sample %in% Brca_norm_plot$sample){
    Brca_norm_plot[which(Brca_norm_plot$sample == sample & Brca_norm_plot$hgnc_symbol == 'WEE1'),]$copy_number_alteration <- WEE1[which(WEE1$Sample.Id == sample),]$WEE1..Putative.copy.number.alterations.from.GISTIC
    
  }
}

table(Brca_norm_plot$copy_number_alteration) 

Brca_norm_plot_wt_cn <- Brca_norm_plot %>% na.omit() 


linear_Regression_gene <- new.env()

Brca_norm_plot_wt_cn[which(Brca_norm_plot_wt_cn$hgnc_symbol == 'CTDP1'),]$hgnc_symbol <- 'FCP1'
Brca_norm_plot_wt_cn[which(Brca_norm_plot_wt_cn$hgnc_symbol == 'MARK3'),]$hgnc_symbol <- 'C-TAK1'
Brca_norm_plot_wt_cn[which(Brca_norm_plot_wt_cn$hgnc_symbol == 'CHEK1'),]$hgnc_symbol <- 'CHK1'


#change the values in numbers 
for (gene in list('FOXM1','CDC25C','CCNB1','CHK1','WEE1','CDK1','C-TAK1','FCP1')){
  
  print(gene)
  plotting <- Brca_norm_plot_wt_cn[which(Brca_norm_plot_wt_cn$hgnc_symbol == gene),]
  
  plotting$copy_number_alteration <- as.factor(plotting$copy_number_alteration)
  
  if (length(levels(plotting$copy_number_alteration)) == 5){
    plotting$copy_number_alteration <- factor(plotting$copy_number_alteration, levels =  c('Deep Deletion','Shallow Deletion', 'Diploid','Gain','Amplification'))
 
    plotting <- transform(plotting,copy_number_alteration = as.integer(copy_number_alteration))
    
    lin_reg_log <- lm(log10(value +1) ~ copy_number_alteration, data = plotting)
    
    f <- (summary(lin_reg_log)$fstatistic)
    pval <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(pval) <- NULL
    
    print(pval)
    if (pval <= 0.01){
    
    linear_Regression_gene[[paste0(gene,'_log')]] <- lin_reg_log
    }

    plotting$copy_number_alteration <- factor(plotting$copy_number_alteration, levels = c(1,2,3,4,5), labels =  c('Deep Deletion','Shallow Deletion', 'Diploid','Gain','Amplification'))
    
    p <- ggplot(plotting, aes(x=factor(copy_number_alteration, levels = c('Deep Deletion','Shallow Deletion', 'Diploid','Gain','Amplification')), y =log10(value +1), fill =copy_number_alteration )) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c('#d7191c',
    '#fdae61','#ffffbf','#abd9e9','#2c7bb6' )) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Copy number status') + ylab('log10(TPM + 1)') + ggtitle(paste0(gene)) +  geom_abline(color = 'darkred',intercept = coef(lin_reg_log)[[1]],slope = coef(lin_reg_log)[[2]])

    print(p) 
  } else if (length(levels(plotting$copy_number_alteration)) == 4 & all( c('Shallow Deletion', 'Diploid','Gain','Amplification') %in% levels(plotting$copy_number_alteration))) {
    plotting$copy_number_alteration <- factor(plotting$copy_number_alteration, levels =  c('Shallow Deletion', 'Diploid','Gain','Amplification'))
    
    plotting <- transform(plotting,copy_number_alteration = as.integer(copy_number_alteration))
    
    lin_reg_log <- lm(log10(value +1) ~ copy_number_alteration, data = plotting)
    
    f <- (summary(lin_reg_log)$fstatistic)
    pval <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(pval) <- NULL
    
    print(pval)
    if (pval <= 0.01){
      
      linear_Regression_gene[[paste0(gene,'_log')]] <- lin_reg_log
    }    
    plotting$copy_number_alteration <- factor(plotting$copy_number_alteration, levels = c(1,2,3,4), labels =  c('Shallow Deletion', 'Diploid','Gain','Amplification'))
    

    p <- ggplot(plotting, aes(x=factor(copy_number_alteration, levels = c('Shallow Deletion', 'Diploid','Gain','Amplification')), y =log10(value +1), fill =copy_number_alteration )) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c('#fdae61','#ffffbf','#abd9e9','#2c7bb6' )) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Copy number status') + ylab('log10(TPM + 1)') + ggtitle(paste0(gene)) +  geom_abline(color = 'darkred',intercept = coef(lin_reg_log)[[1]],slope = coef(lin_reg_log)[[2]])

    print(p)
  } else if(length(levels(plotting$copy_number_alteration)) == 4 & all( c('Deep Deletion','Shallow Deletion', 'Diploid','Gain') %in% levels(plotting$copy_number_alteration))){
    plotting$copy_number_alteration <- factor(plotting$copy_number_alteration, levels =  c('Deep Deletion','Shallow Deletion', 'Diploid','Gain'))
    
    plotting <- transform(plotting,copy_number_alteration = as.integer(copy_number_alteration))
    
    lin_reg_log <- lm(log10(value +1) ~ copy_number_alteration, data = plotting)
    
    f <- (summary(lin_reg_log)$fstatistic)
    pval <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(pval) <- NULL
    
    print(pval)
    if (pval <= 0.01){
      
      linear_Regression_gene[[paste0(gene,'_log')]] <- lin_reg_log
    }
    plotting$copy_number_alteration <- factor(plotting$copy_number_alteration, levels = c(1,2,3,4), labels = c('Deep Deletion','Shallow Deletion', 'Diploid','Gain'))
  
    p <- ggplot(plotting, aes(x=factor(copy_number_alteration, levels = c('Deep Deletion','Shallow Deletion', 'Diploid','Gain')), y =log10(value +1), fill =copy_number_alteration )) + geom_boxplot() +scale_fill_manual(values = c('#d7191c','#fdae61','#ffffbf','#abd9e9')) +  theme_minimal() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('Copy number status') + ylab('log10(TPM + 1)') + ggtitle(paste0(gene)) +  geom_abline(color = 'darkred',intercept = coef(lin_reg_log)[[1]],slope = coef(lin_reg_log)[[2]])
 
    print(p)
    
    }
  
}

#filtered already for p.value <
CDK1_delta <- predict.lm(linear_Regression_gene[["CDK1_log"]],data.frame(copy_number_alteration = c(1,2,3,4,5)))
CHK1_delta <- predict.lm(linear_Regression_gene[["CHK1_log"]],data.frame(copy_number_alteration = c(1,2,3,4)))
FCP1_delta <- predict.lm(linear_Regression_gene[["FCP1_log"]],data.frame(copy_number_alteration = c(1,2,3,4,5)))
FOXM1_delta <- predict.lm(linear_Regression_gene[["FOXM1_log"]],data.frame(copy_number_alteration = c(1,2,3,4)))
CTAK1_delta <- predict.lm(linear_Regression_gene[["C-TAK1_log"]],data.frame(copy_number_alteration = c(1,2,3,4,5)))
WEE1_delta <- predict.lm(linear_Regression_gene[["WEE1_log"]],data.frame(copy_number_alteration = c(1,2,3,4,5)))

# CCNB1 and cdc25c doesn't seems to follow a linear distribution, so not present 

CDK1_diffrences <- data.frame(Gene = 'CDK1',Copy_number = c('Deep Deletion','Shallow Deletion', 'Diploid','Gain','Amplification'), values = NA, linear_reg_values =CDK1_delta )
CHK1_diffrences <- data.frame(Gene = 'CHK1',Copy_number = c('Deep Deletion','Shallow Deletion', 'Diploid','Gain'), values = NA, linear_reg_values =CHK1_delta )
FCP1_diffrences <- data.frame(Gene = 'FCP1',Copy_number =c('Deep Deletion','Shallow Deletion', 'Diploid','Gain','Amplification'), values = NA, linear_reg_values =FCP1_delta )
FOXM1_diffrences <- data.frame(Gene = 'FOXM1',Copy_number =c('Shallow Deletion', 'Diploid','Gain','Amplification'), values = NA, linear_reg_values =FOXM1_delta )
CTAK1_diffrences <- data.frame(Gene = 'C-TAK1',Copy_number =c('Deep Deletion','Shallow Deletion', 'Diploid','Gain','Amplification'), values = NA, linear_reg_values =CTAK1_delta )
WEE1_diffrences <- data.frame(Gene = 'WEE1',Copy_number =c('Deep Deletion','Shallow Deletion', 'Diploid','Gain','Amplification'), values = NA, linear_reg_values =WEE1_delta )


ratios_genes_copynumb <- do.call("rbind", list(WEE1_diffrences, CTAK1_diffrences, FOXM1_diffrences,FCP1_diffrences,CHK1_diffrences,CDK1_diffrences))

# insert ratios 
for (gene_copy in ratios_genes_copynumb$Gene){
  momentary <- ratios_genes_copynumb[which(ratios_genes_copynumb$Gene == gene_copy),]
  for (status in momentary$Copy_number){
    if(status == 'Diploid'){
      ratios_genes_copynumb[which(ratios_genes_copynumb$Gene == gene_copy & ratios_genes_copynumb$Copy_number == status),]$values <- 500
    }
    else{
      ratios_genes_copynumb[which(ratios_genes_copynumb$Gene == gene_copy & ratios_genes_copynumb$Copy_number == status),]$values <- (500/momentary[which(momentary$Copy_number== 'Diploid'),]$linear_reg_values) * momentary[which(momentary$Copy_number== status),]$linear_reg_values
    }
  }
}


write_csv(ratios_genes_copynumb,file = '../Files_generated/final_rates_copynumbers.csv')

table_patient_specific <- subset(Brca_norm_plot_wt_cn,select=c('sample','hgnc_symbol','copy_number_alteration'))

for(names in unique(table_patient_specific$hgnc_symbol)){
  if(names == 'CCNB1' || names == 'CDC25C'){
   table_patient_specific[which(table_patient_specific$hgnc_symbol == names),]$copy_number_alteration <- NA
  }
}


table_patient_specific <- table_patient_specific %>%  na.omit()

# check if some samples miss genes
count_patient = 0
right_count = 0
wrong_count = 0
for(patient in unique(table_patient_specific$sample)){
  selcted_ones <- table_patient_specific[which(table_patient_specific$sample == patient),]
  genes_retrived <- selcted_ones$hgnc_symbol
  count_patient = count_patient + 1
  if(length(genes_retrived) == 7){
    right_count = right_count + 1
  }else{
    wrong_count = wrong_count + 1
  }
}



write_csv(table_patient_specific,'../Files_generated/table_patient_specific.csv')  


