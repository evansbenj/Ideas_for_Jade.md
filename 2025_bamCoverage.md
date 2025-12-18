# bamCoverage 

directory:
```
/home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/deepTools/deeptools
```
bamCoverage is a tool in the package "deepTools" that calculates coverage in windows from bam files. Here is a script that runs it:
```
/home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/ben_scripts/2025_bamCoverage.sh
```

```sh
#!/bin/sh
#SBATCH --job-name=bamCoverage
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=bamCoverage.%J.out
#SBATCH --error=bamCoverage.%J.err
#SBATCH --account=rrg-ben

module load StdEnv/2023 python/3.13.2

/home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/deepTools/deeptools/vcf_env/bin/bamCoverage -b ${1} --normalizeUsing RPKM --outFileFormat bedgraph --binSize 100000 --ignoreDuplicates --minMappingQuality 30 -o ${1}_bamCoverage.bw
```
# Plot
```
library (ggplot2)
library(tidyverse)
library(reshape2) # this facilitates the overlay plot
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2024_cliv_allo_WGS/2025_allo_WGS_depth_in_windowz")
# read in the data 
options(scipen=999)
# these coverages are counts from bamCoverage. I am pretty sure that they are the number of reads
# that map to each window (in this case each is 100,000 bp long)
# to convert to coverage per position, multiply the number of reads by the average length of
# a pair of reads (~300bp) and then divide by the window size (100,000)
# this provides ~22X coverage for an average count of 7500 per window, which seems correct

# declare an empty df
#standardized_depth <- data.frame(matrix(ncol = 10, nrow = 0))
#colnames(standardized_depth) <- c('BJE3485', 'BJE3486', 'BJE3487','BJE3488', 'BJE3495', 'BJE3496','BJE3501', 'BJE3502', 'Z23732','Z23738')

depth_per_site_BJE3485 <- read.table("BJE3485_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_BJE3485) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_BJE3485)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_BJE3485$depth));x
stdev <- sd(depth_per_site_BJE3485$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_BJE3485$depth,
                                            function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_BJE3485$standardized_depth <- sapply(depth_per_site_BJE3485$depth,
                                                         function(r) ifelse ((is.numeric(r)),r/x,NA))


depth_per_site_BJE3486 <- read.table("BJE3486_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_BJE3486) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_BJE3486)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_BJE3486$depth));x
stdev <- sd(depth_per_site_BJE3486$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_BJE3486$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_BJE3486$standardized_depth <- sapply(depth_per_site_BJE3486$depth,
                                                    function(r) ifelse ((is.numeric(r)),r/x,NA))

depth_per_site_BJE3487 <- read.table("BJE3487_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_BJE3487) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_BJE3487)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_BJE3487$depth));x
stdev <- sd(depth_per_site_BJE3487$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_BJE3487$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_BJE3487$standardized_depth <- sapply(depth_per_site_BJE3487$depth,
                                                    function(r) ifelse ((is.numeric(r)),r/x,NA))


depth_per_site_BJE3488 <- read.table("BJE3488_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_BJE3488) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_BJE3488)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_BJE3488$depth));x
stdev <- sd(depth_per_site_BJE3488$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_BJE3488$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_BJE3488$standardized_depth <- sapply(depth_per_site_BJE3488$depth,
                                                    function(r) ifelse ((is.numeric(r)),r/x,NA))


depth_per_site_BJE3495 <- read.table("BJE3495_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_BJE3495) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_BJE3495)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_BJE3495$depth));x
stdev <- sd(depth_per_site_BJE3495$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_BJE3495$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_BJE3495$standardized_depth <- sapply(depth_per_site_BJE3495$depth,
                                                    function(r) ifelse ((is.numeric(r)),r/x,NA))

depth_per_site_BJE3496 <- read.table("BJE3496_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_BJE3496) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_BJE3496)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_BJE3496$depth));x
stdev <- sd(depth_per_site_BJE3496$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_BJE3496$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_BJE3496$standardized_depth <- sapply(depth_per_site_BJE3496$depth,
                                                    function(r) ifelse ((is.numeric(r)),r/x,NA))

depth_per_site_BJE3501 <- read.table("BJE3501_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_BJE3501) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_BJE3501)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_BJE3501$depth));x
stdev <- sd(depth_per_site_BJE3501$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_BJE3501$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_BJE3501$standardized_depth <- sapply(depth_per_site_BJE3501$depth,
                                                    function(r) ifelse ((is.numeric(r)),r/x,NA))

depth_per_site_BJE3502 <- read.table("BJE3502_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_BJE3502) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_BJE3502)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_BJE3502$depth));x
stdev <- sd(depth_per_site_BJE3502$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_BJE3502$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_BJE3502$standardized_depth <- sapply(depth_per_site_BJE3502$depth,
                                                    function(r) ifelse ((is.numeric(r)),r/x,NA))



depth_per_site_Z23732 <- read.table("Z23732_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_Z23732) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_Z23732)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_Z23732$depth));x
stdev <- sd(depth_per_site_Z23732$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_Z23732$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_Z23732$standardized_depth <- sapply(depth_per_site_Z23732$depth,
                                                   function(r) ifelse ((is.numeric(r)),r/x,NA))


depth_per_site_Z23738 <- read.table("Z23738_sorted_rg_dedup.bam_bamCoverage.bw", header = F)
colnames(depth_per_site_Z23738) <- c("contig","start","stop","depth")
#table(as.character(sapply(depth_per_site, class)))
dim(depth_per_site_Z23738)


# Standardize the depth by the mean depth per sample
# first calculate the mean count with outliers included
x <- mean(na.omit(depth_per_site_Z23738$depth));x
stdev <- sd(depth_per_site_Z23738$depth, na.rm = TRUE); stdev

# now use this mean to get rid of outliers
depth_per_sitez <- sapply(depth_per_site_Z23738$depth,
                          function(r) ifelse( (is.numeric(r))&&(r<(stdev*3)),r,NA))
# recalculate mean
x <- mean(na.omit(depth_per_sitez));x
# standardize
depth_per_site_Z23738$standardized_depth <- sapply(depth_per_site_Z23738$depth,
                                                   function(r) ifelse ((is.numeric(r)),r/x,NA))


library(dplyr)
# malez_first
df3 <- full_join(depth_per_site_BJE3485, depth_per_site_BJE3486, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(df3) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth')
df4 <- full_join(df3, depth_per_site_BJE3495, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(df4) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth',
                   'BJE3495_depth','BJE3495_stdepth')
df5 <- full_join(df4, depth_per_site_BJE3496, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(df5) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth',
                   'BJE3495_depth','BJE3495_stdepth','BJE3496_depth','BJE3496_stdepth')
df6 <- full_join(df5, depth_per_site_Z23738, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(df6) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth',
                   'BJE3495_depth','BJE3495_stdepth','BJE3496_depth','BJE3496_stdepth',
                   'Z23738_depth','Z23738_stdepth')

# now add femalez
df7 <- full_join(df6, depth_per_site_BJE3487, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(df7) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth',
                   'BJE3495_depth','BJE3495_stdepth','BJE3496_depth','BJE3496_stdepth',
                   'Z23738_depth','Z23738_stdepth','BJE3487_depth','BJE3487_stdepth')

df8 <- full_join(df7, depth_per_site_BJE3488, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(df8) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth',
                   'BJE3495_depth','BJE3495_stdepth','BJE3496_depth','BJE3496_stdepth',
                   'Z23738_depth','Z23738_stdepth','BJE3487_depth','BJE3487_stdepth',
                   'BJE3488_depth','BJE3488_stdepth')

df9 <- full_join(df8, depth_per_site_BJE3501, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(df9) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth',
                   'BJE3495_depth','BJE3495_stdepth','BJE3496_depth','BJE3496_stdepth',
                   'Z23738_depth','Z23738_stdepth','BJE3487_depth','BJE3487_stdepth',
                   'BJE3488_depth','BJE3488_stdepth','BJE3501_depth','BJE3501_stdepth')

df10 <- full_join(df9, depth_per_site_BJE3502, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(df10) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth',
                   'BJE3495_depth','BJE3495_stdepth','BJE3496_depth','BJE3496_stdepth',
                   'Z23738_depth','Z23738_stdepth','BJE3487_depth','BJE3487_stdepth',
                   'BJE3488_depth','BJE3488_stdepth','BJE3501_depth','BJE3501_stdepth',
                   'BJE3502_depth','BJE3502_stdepth')

all_of_em <- full_join(df10, depth_per_site_Z23732, by=c('contig'='contig', 'start'='start', 'stop'='stop'))
colnames(all_of_em) <- c('contig','start','stop','BJE3485_depth','BJE3485_stdepth','BJE3486_depth','BJE3486_stdepth',
                    'BJE3495_depth','BJE3495_stdepth','BJE3496_depth','BJE3496_stdepth',
                    'Z23738_depth','Z23738_stdepth','BJE3487_depth','BJE3487_stdepth',
                    'BJE3488_depth','BJE3488_stdepth','BJE3501_depth','BJE3501_stdepth',
                    'BJE3502_depth','BJE3502_stdepth','Z23732_depth','Z23732_stdepth')

# calculate the F-M difference in mean depth 
all_of_em$mean_diff_M_minus_F <- rowMeans(all_of_em[ , c(5,7,9,11,13)], na.rm=TRUE) -
  rowMeans(all_of_em[ , c(15,17,19,21,23)], na.rm=TRUE)
all_of_em$mean_malez <- rowMeans(all_of_em[ , c(5,7,9,11,13)], na.rm=TRUE) 
all_of_em$mean_femalez <- rowMeans(all_of_em[ , c(15,17,19,21,23)], na.rm=TRUE)

all_of_em_noNAs <- all_of_em[complete.cases(all_of_em), ]


#all_of_em_noNAs %>% select(contig, mean_diff_F_minus_M_notads, start) %>%
#  pivot_longer(., cols = c(mean_diff_F_minus_M_notads), names_to = "Var", values_to = "Val") %>%
#  ggplot(aes(x = Var, y = Val, fill=contig)) +
#  geom_point()


# Manhattan plot function ----
library(tidyverse)
#install.packages("qqman")
library(qqman)
library(dplyr)



### BELOW works using lattice graphics instead of ggplot
# but there are issues with adding rectangles, etc
# but otherwise it looks pretty good
library(lattice)
library(plyr)
library(latticeExtra)
# manhattan plot function ----
# This script is from:
# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab="Mean Depth Diff (M-F)",
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(pvalue,thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- pvalue
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=45,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0)))+.1;
    miny<-floor(min(y, ifelse(!is.na(sig.level), sig.level, 0)))-.1;
    A$ylim=c(miny,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=sig.level, lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}
# chromosome plot function ----

chromosome.plot<-function(chr, pos, pvalue, 
                          sig.level=NA, annotate=NULL, ann.default=list(),
                          should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                          xlab="Mb", ylab="Mean depth diff (M-F)",
                          col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fonte=18, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(pvalue,thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- pvalue
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 #at=((posmax+posmin)/2+posshift),
                 labels=T, 
                 ticks=T, rot=0,
                 text.cex=1,
                 #tck = seq(0,200,50),
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    } else if(side=="left") {
      panel.axis(side=side, outside=T,
                 #at=((posmax+posmin)/2+posshift),
                 labels=T, 
                 ticks=T, rot=0,
                 text.cex=1,
                 #tck = seq(0,200,50),
                 check.overlap=F
      )
    } else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0)))+.5;
    miny<-floor(min(y, ifelse(!is.na(sig.level), sig.level, 0)))-.1;
    A$ylim=c(miny,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=sig.level, lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}

theme.novpadding <- list(
  layout.heights = list(
    top.padding = 0,
    main.key.padding = 0,
    key.axis.padding = 0,
    axis.xlab.padding = 0,
    xlab.key.padding = 0,
    key.sub.padding = 0,
    bottom.padding = 0
  ),
  layout.widths = list(
    left.padding = 0,
    key.ylab.padding = 0,
    ylab.axis.padding = 0,
    axis.key.padding = 0,
    right.padding = 0
  )
)

# chromosome plot no thin function ----

chromosome.plot.no.thin<-function(chr, pos, pvalue, 
                                  sig.level=NA, annotate=NULL, ann.default=list(),
                                  should.thin=F, thin.pos.places=2, thin.logp.places=2, 
                                  xlab="Mb", ylab="Mean Depth Diff (M-F)",
                                  col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fonte=18, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(pvalue,thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- pvalue
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 #at=((posmax+posmin)/2+posshift),
                 labels=T, 
                 ticks=T, rot=0,
                 text.cex=1,
                 #tck = seq(0,200,50),
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    } else if(side=="left") {
      panel.axis(side=side, outside=T,
                 #at=((posmax+posmin)/2+posshift),
                 labels=T, 
                 ticks=T, rot=0,
                 text.cex=1,
                 #tck = seq(0,200,50),
                 check.overlap=F
      )
    } else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0)))+.5;
    miny<-floor(min(y, ifelse(!is.na(sig.level), sig.level, 0)))-.1;
    A$ylim=c(miny,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=sig.level, lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}

# get rid of tig00007285 because this is mtDNA
all_of_em_noNAs <- all_of_em_noNAs[all_of_em_noNAs$contig != "tig00007285",]
# change contig names
all_of_em_noNAs <- as.data.frame(sapply(all_of_em_noNAs,function(x) {x <- gsub("tig0000","c",x)}))
all_of_em_noNAs <- as.data.frame(sapply(all_of_em_noNAs,function(x) {x <- gsub("tig000","c",x)}))
# get variables in proper format
all_of_em_noNAs$start <- as.numeric(all_of_em_noNAs$start)
all_of_em_noNAs$mean_diff_M_minus_F <- as.numeric(all_of_em_noNAs$mean_diff_M_minus_F)
all_of_em_noNAs$mean_femalez <- as.numeric(all_of_em_noNAs$mean_femalez)
all_of_em_noNAs$mean_malez <- as.numeric(all_of_em_noNAs$mean_malez)




# plot ----
M_minus_F_allo_depth_plot <- manhattan.plot(factor(all_of_em_noNAs$contig), all_of_em_noNAs$start, 
                                            all_of_em_noNAs$mean_diff_M_minus_F)#,
                                          # ylim= c(0,12),xlab = "",
                                         # par.settings = theme.novpadding)#,
                                          #panel=function(x, y, ...){
                                            # panel.rect(xleft=695, ybottom=0,
                                            #           xright=733, ytop=10, alpha=0.3, col="light blue")
                                          #panel.text(275,9,labels=expression(italic("X. allofraseri")),fontsize=14)})

jpeg("./M_minus_F_allo_depth_plot.jpg",w=500, h=3.0, units ="in", bg="transparent", res = 200)
  M_minus_F_allo_depth_plot
dev.off()



SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c11059",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c10978",] # carries foxi2.L
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c10680",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c10638",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c10588",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c10461",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c10385",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c10342",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c10341",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c4422",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c7519",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c7399",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c9267",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c8579",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c7537",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c6777",]
SexChr <- all_of_em_noNAs[all_of_em_noNAs$contig == "c1724",]




# narrow down region
#SexChr <- SexChr[((SexChr$pos > 2300000)&(SexChr$pos < 2700000)),]
Contig_focus_plot <- chromosome.plot.no.thin(factor(SexChr$contig),SexChr$start,SexChr$mean_diff_M_minus_F)
                                              

jpeg("./Contig_focus_plot_c7399_focus.jpg",w=7, h=3.0, units ="in", bg="transparent", res = 200)
  Contig_focus_plot
dev.off()

View(SexChr)
```
