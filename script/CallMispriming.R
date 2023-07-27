## Identify mispriming ITRseq peaks

# fbam <- "/project/gtplab/data_analysis/C598-ITR-run4-5-rerun/22-130_Heart_Ventricle_R_ref_map_filtered_sorted.bam"
# fbed <- "/project/gtplab/data_analysis/C598-ITR-run4-5-rerun/22-130_Heart_Ventricle_R_ref_peak.bed"
# fbam <- "/project/gtplab/data_analysis/C592-ITR/22-130_Liver_Caudate_ref_map_filtered_sorted_dedup_novec.bam"
# fbed <- "/project/gtplab/data_analysis/C592-ITR/22-130_Liver_Caudate_ref_peak.bed"
# fbam <- "/project/gtplab/data_analysis/C592-ITR/22-173_Liver_Biopsy_2020-007-3029-001_ref_map_filtered_sorted_dedup_novec.bam"
# fbed <- "/project/gtplab/data_analysis/C592-ITR/22-173_Liver_Biopsy_2020-007-3029-001_ref_peak.bed"
# primer.seq <- "ACAAGGAACCCCTAGTGATGGAGTTGGCC";
# seed.seq <- "ACTCCCTC";
# max.mismatch <- 0;
# min.total=1;
# min.strand=0;

CallMispriming <- function(fbam, fbed, seed.seq, max.mismatch=0, min.total=1, min.strand=0) {
  # fbam          Name of bam file with aligned reads
  # fbed          Name of bed file with genomic location of the ITRseq peaks
  # seed.seq      Seed sequence in the ITR to be searched in reads
  # max.mismatch  Max number of mismatches allowed when match seed sequence to reads
  # min.total     Min of total selected reads not to call a peak mispriming
  # min.strand    Min of selected reads on both strands not to call a peak mispriming
  
  require(GenomicAlignments);
  
  ## Also use reverse-complement sequence as seed
  seed.seq <- DNAStringSet(seed.seq);
  rev.seq <- reverseComplement(seed.seq);
  sseq <- c(seed.seq, rev.seq);
  
  if (file.info(fbed[1])$size == 0) NULL else {
    
    bed <- read.csv(fbed, header = FALSE, sep='\t', stringsAsFactors = FALSE);
    gr  <- GRanges(bed[, 1], IRanges(bed[, 2], bed[, 3]));
    flg <- scanBamFlag(isPaired = TRUE, isSecondaryAlignment = FALSE);
    
    print(fbed);
    
    cnt <- lapply(1:length(gr), function(i) { print(i);
      aln <- readGAlignmentPairs(fbam, param = ScanBamParam(what=c('qname', 'cigar', 'seq'), which=gr[i], flag=flg));
      
      read2 <- aln[countOverlaps(aln, gr[i], maxgap = 1)>0]@last;
      
      if (length(read2) == 0) rep(0, 6) else {
        
        names(read2) <- 1:length(read2);
        
        ## Match seed sequence 
        seed.match <- lapply(1:length(sseq), function(j) {
          s <- sseq[[j]];
          m <- vmatchPattern(s, read2@elementMetadata$seq, max.mismatch = max.mismatch);
          q <- rep(names(read2), elementNROWS(m));
          m <- unlist(m);
          GRanges(q, m);
        });
        seed.match <- suppressWarnings(do.call('c', seed.match));
        
        ## Region within read matching genomic sequence according to CIGAR
        genome.match <- cigarRangesAlongQuerySpace(cigar(read2), ops = 'M');
        genome.match <- GRanges(rep(names(read2), elementNROWS(genome.match)), unlist(genome.match));
        
        match <- seed.match[countOverlaps(seed.match, genome.match, type='within')==0];
        sel  <- read2[unique(as.vector(seqnames(match)))];
        
        str0 <- as.vector(strand(read2));
        str1 <- as.vector(strand(sel));
        c(length(read2), length(read2[str0=='+']), length(read2[str0=='-']), length(sel), length(sel[str1=='+']), length(sel[str1=='-']));
      }
    });
    cnt <- do.call('rbind', cnt);
    colnames(cnt) <- c('Total', 'Plus', 'Minus', 'Seeded', 'Seeded_Plus', 'Seeded_Minus');
    rownames(cnt) <- paste0(bed[, 1], ':', bed[, 2], '-', bed[, 3]);

    score <- bed[,5];
    
    out <- cbind(cnt, Score=score, Mispriming=1-as.integer(cnt[, 4]>=min.total & cnt[, 5]>=min.strand & cnt[, 6]>=min.strand));
    out[rowSums(cnt)==0, 7] <- 0;
    
    invisible(out);
  }
}
