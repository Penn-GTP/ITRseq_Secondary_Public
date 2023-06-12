library(stringr)
library(BSgenome, lib.loc = "/appl/R-4.2/lib64/R/library")
library(BSgenome.Mmulatta.UCSC.rheMac10)

file_C592 <- read_csv("file_C592.csv")

########################
# Mispriming

seed.seq <- "ACTCCCTC"

mispriming <- mapply(CallMispriming, fbam = file_C592$bam, fbed = file_C592$bed, seed.seq = seed.seq)
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

anno <- mapply(ReformatPeakAnno, file_C592$anno)
anno_list <- list()
smp_names <- list()
for (i in 1:length(anno)){
  temp <- data.frame(anno[i])
  tryCatch(
    expr={
  colnames(temp) <- c("Peak", "Chrom", "Start", "End", "Read", "Score", "Gene_ID", "Gene_Name", "Gene_Biotype", "Gene_Description", "Gene_TSS", "Gene_TSS2", "Strand", "Overlap", "Dist2TSS", "CDS", "Promotor", "5UTR", "3UTR", "Intragenic", "Intergenic", "Representative")
    },
  error = function(e){})
  
  anno_list[[i]] <- temp
  smp_names[[i]] <- str_sub(unlist(strsplit(names(anno[i]), '/'))[6], 1, -19)
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
