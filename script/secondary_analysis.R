library(stringr)
library(BSgenome, lib.loc = "/appl/R-4.2/lib64/R/library")
library(BSgenome.Mmulatta.UCSC.rheMac10)

wd <- unlist(str_split(getwd(), '/'))
parent_dir <- paste0(wd[2:length(wd)-1], collapse = '/')

bams <- list.files(path = '../', pattern = '.*_ref_map_filtered_sorted_dedup_novec.bam$')
abspath_bams <- paste(parent_dir, bams, sep = '/')

beds <- list.files(path = '../', pattern = '.*_ref_peak.bed$')
abspath_beds <- paste(parent_dir, beds, sep = '/') 

annos <- list.files(path = '../', pattern = '.*_ref_peak_anno.tsv$')
abspath_annos <- paste(parent_dir, annos, sep = '/') 

name_idx <- mapply(lcprefix, bams, beds) - 5
ids <- substr(bams, 1, name_idx)

file_paths <- data.frame(bam = abspath_bams, bed = abspath_beds, anno = abspath_annos, row.names = ids)
write.csv(file_paths, 'filepaths.csv')

#file_C592 <- read_csv("file_C592.csv")

### Remove On-Target ###

library(bedtoolsr)

off_target_beds <- mapply(bt.window, a=file_paths$bed, b='ARCUS_target_Mmul10.bed', w=500, v=TRUE, SIMPLIFY = FALSE)

for (i in 1:length(off_target_beds)){
  write_tsv(data.frame(off_target_beds[i]), paste0('offtarget_beds/offtarget_', unlist(strsplit(names(off_target_beds[i]), '/'))[6]), col_names = FALSE)
}

########################
# Mispriming


seed.seq <- "ACTCCCTC"

mispriming <- mapply(CallMispriming, fbam = file_paths$bam, fbed = paste0( 'offtarget_beds/' ,list.files(path = 'offtarget_beds', pattern = '.*bed$')), seed.seq = seed.seq)
mispriming_list <- list()
for (i in 1:length(mispriming)){
  temp <- data.frame(mispriming[i])
  tryCatch(
    expr = {
      colnames(temp) <- c('Total','Plus','Minus','Seeded','Plus_seeded','Minus_seeded','Mispriming')
    },
    error = function(e){ #do nothing
    })
  
  mispriming_list[[i]] <- temp
  write.csv(temp, paste('mispriming/',str_sub(unlist(strsplit(names(mispriming[i]), '/'))[6], 1, -41),'_mispriming.csv', sep=''))
}

##############################
# Annotation

anno <- mapply(ReformatPeakAnno, file_paths$anno)
anno_list <- list()
smp_names <- list()
for (i in 1:length(anno)){
  temp <- data.frame(anno[i])
  tryCatch(
    expr={
      colnames(temp) <- c("Peak", "Chrom", "Start", "End", "Read", "Score", "Gene_ID", "Gene_Name", "Gene_Biotype", "Gene_Description", "Gene_TSS", "Gene_TTS", "Strand", "Overlap", "Dist2TSS", "CDS", "Promoter", "5UTR", "3UTR", "Intragenic", "Intergenic", "Representative")
    },
    error = function(e){})
  
  
  smp_names[[i]] <- str_sub(unlist(strsplit(names(anno[i]), '/'))[6], 1, -19)
  
  tryCatch(
    expr={
      temp$sample <- smp_names[[i]]    },
    error = function(e){})
  
  anno_list[[i]] <- temp
  
  write.csv(temp, paste('annotation/',str_sub(unlist(strsplit(names(anno[i]), '/'))[6], 1, -19),'_anno.csv', sep=''))
}
#anno <- lapply(anno , function(x){ row.names(x)<-as.character(x$Peak); x}) #set the peak to rownames

###################
# Match to PWM

pwm <- as.matrix(readRDS('PWM.rds'))
for (i in 1:length(mispriming_list)){
  print(i)
  mispriming_list[[i]][['Mispriming']] = as.logical(mispriming_list[[i]][['Mispriming']]) #changes misprimed peaks from 1 to TRUE
  masked_anno <- anno_list[[i]][!mispriming_list[[i]][['Mispriming']],] # only keeps FALSE in mispriming
  tryCatch(
    exp = {
      match <- MatchGenomeToPWM(Mmulatta, paste0('chr', masked_anno$Chrom), masked_anno$Start, masked_anno$End, pwm)
      masked_anno$PWM_Score <- match$Score
      write.csv(masked_anno, paste0('PWM_match/', smp_names[[i]], '.csv'))
    },
    error = function(e){ 
      print(paste(smp_names[[i]], 'has no non-mispriming peaks'))
    }
  )
}

####################
# Concat result dataframes

library(dplyr)
library(tidyverse)
library(bedr)

score_dfs_list <- list.files(path = 'PWM_match')


score_dfs<-mapply(read.csv, file = paste('PWM_match', score_dfs_list, sep = '/'),  SIMPLIFY = FALSE)

score_dfs_fixed <- lapply(score_dfs, function(df) {
  transform(df, Chrom = as.character(Chrom))})
combined <- bind_rows(score_dfs_fixed)

tempbed <- data.frame(chr=combined$Chrom, start=combined$Start, end=combined$End, name=combined$sample, score=combined$Score, strand=combined$Strand)


out_name <- '2020-007-nonliver-no-mispriming-peaks.bed'
write_tsv(tempbed, out_name, col_names = FALSE)
out.sort <- bt.sort(out_name)
out.merged <- bt.merge(out.sort, d=44, c='4,5,6', o='collapse,sum,distinct')

write_tsv(out.merged, out_name, col_names = FALSE)

