## Find the best match to a given PWM in given location of a reference genome and calculate its matching score
## Best match to PWM in a sequence

################################################################################################
################################################################################################
CalcPWMScore <- function(pwm, seq, nbase=1) {
  smx <- sum(apply(pwm, 2, max));
  smn <- sum(apply(pwm, 2, min));
  
  seq <- as.character(seq);
  
  bas <- sapply(1:(nchar(seq)-nbase+1), function(i) substr(seq, i, i+nbase-1));
  
  mtx <- pwm[bas, ];
  sum <- sapply(1:(nrow(mtx)-ncol(mtx)+1), function(i) sum(diag(mtx[i:(i+ncol(mtx)-1), ])));
  
  (sum-smn)/(smx-smn);
}

################################################################################################
################################################################################################
BestPwmMatch <- function(s, pwm, nbase=1) {
  require(Biostrings);
  
  s   <- DNAString(s);
  v   <- reverseComplement(s);
  sc1 <- CalcPWMScore(pwm, s, nbase=nbase);
  sc2 <- CalcPWMScore(pwm, v, nbase=nbase);
  # sc1 <- sapply(1:(s@length-ncol(pwm)+1), function(i) PWMscoreStartingAt(pwm, s, i));
  # sc2 <- sapply(1:(s@length-ncol(pwm)+1), function(i) PWMscoreStartingAt(pwm, v, i));
  sc  <- c(sc1, sc2);
  str <- rep(c(1, -1), c(length(sc1), length(sc2)));
  pos <- c(1:length(sc1), s@length - (1:length(sc2)) - (ncol(pwm)-1) - (nbase - 1) + 1);
  dst <- abs(rep(1:length(sc1), 2) - (s@length/2 - ncol(pwm)/2));
  whh <- which(sc == max(sc));
  if (length(whh) > 1) {
    d <- dst[whh];
    whh <- sample(whh[d==min(d)], 1);
  }
  c(position=pos[whh], strand=str[whh], score=sc[whh]) 
}

MatchGenomeToPWM <- function(gref, chrom, start, end, pwm, ext=ncol(pwm)) {
  # gref    BSgenome object of the reference genome
  # chrom, start, end   Genomic location in the reference genome
  # pwm     PWM of a motif to be matched; must be count of ACGT bases at each position
  # ext     Number of extra bases to be added to the genomic location; on both ends
  
  require(BSgenome);
  require(GenomicRanges);
  
  ###################################################################################################
  
  summarizeBestMatch <- function(bst0, pwm, site0, nbase=1) {
    match <- t(sapply(as.character(bst0$full), function(s) bestPwmMatch(s, pwm, nbase=nbase)));
    
    bst1 <- bst0;
    bst1@ranges <- IRanges(bst1$peak_start + match[, 1] - 1, bst1$peak_start + match[, 1] + width(bst0) - 2);
    bst1@strand <- Rle(factor(c('-', '+')[1 + pmax(0, match[, 2])], levels=c('+', '-', '*')));
    bst1$score  <- match[, 3];
    
    sfll <- bst1$full;
    ssub <- substr(sfll, start(bst1)-bst1$peak_start+1, end(bst1)-bst1$peak_start+1)
    ssub <- DNAStringSet(ssub);
    ssub[as.character(strand(bst1))=='-'] <- reverseComplement(ssub[as.character(strand(bst1))=='-']);
    
    bst1$seq    <- DNAStringSet(ssub);
    names(bst1) <- names(bst0);
    
    
    bas0 <- strsplit(oseq, '')[[1]];
    bas1 <- do.call('rbind', strsplit(as.character(bst1$seq), ''));
    bst1$nmatch <- apply(bas1, 1, function(b) length(b[b==bas0]));
    
    dst1 <- rep(0, length(bst1));
    dst1[end(bst1) < start(site0)] <- (start(site0) - end(bst1))[end(bst1) < start(site0)];
    dst1[start(bst1) > end(site0)] <- (start(bst1) - end(site0))[start(bst1) > end(site0)];
    bst1$distance <- dst1;
    
    bst1;
  }
  ###################################################################################################
  
  pwm <- PWM(pwm);
  ext <- max(0, ext[1]);
  gr  <- GRanges(chrom, IRanges(start-ext, end+ext));
  seq <- getSeq(gref, gr);
  hit <- t(sapply(seq, function(s) BestPwmMatch(s, pwm, 1)));
  
  tbl <- data.frame(Chrom=chrom, Start=start(gr)+hit[, 1]-1, End=start(gr)+hit[, 1]-1 + (ncol(pwm)-1), 
                    Strand=hit[, 2], Score=hit[, 3], stringsAsFactors = FALSE);
  sub <- substr(seq, hit[, 1], hit[, 1]+ncol(pwm)-1);
  if (nrow(tbl[tbl$Strand==-1, , drop=FALSE]) > 0) 
    sub[tbl$Strand==-1] <- as.character(reverseComplement(DNAStringSet(sub[tbl$Strand==-1])));
  tbl$Seq <- sub;
  
  invisible(tbl);
  
}
