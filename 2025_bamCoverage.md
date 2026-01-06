# bamCoverage 

directory:
```
/home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/deepTools/deeptools
```
bamCoverage is a tool in the package "deepTools" that calculates coverage in windows from bam files. Here is a script that runs it:
```
/home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/ben_scripts/2025_bamCoverage.sh
```
on bam files in this directory:
```
/home/ben/projects/rrg-ben/ben/2024_cliv_allo_WGS/2025_more_allo_WGS/
```
and this one:
```
/home/ben/projects/rrg-ben/ben/2024_cliv_allo_WGS/fq/2024_allo
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
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2024_cliv_allo_WGS/2025_allo_WGS_depth_in_windowz/depth_with_bubbles_normalized")
# read in the data 
rm(list=ls()) # removes all variables

options(scipen=999)

sample_vector <- c("female_allo_17A7","female_allo_17AB","female_allo_17C2","female_allo_17E2","female_BJE3487","female_BJE3488","female_BJE3501","female_BJE3502","female_camF1","female_Z23698","female_Z23702","female_Z23721","female_Z23726","female_Z23732","female_Z23733","male_allo_17A1","male_allo_1873","male_allo_1876","male_allo_18A1","male_allo_18AF","male_allo_190F","male_BJE3485","male_BJE3486","male_BJE3495","male_BJE3496","male_camM5","male_Z23701","male_Z23738","male_Z23739")
#sample_vector_onlymainland <- c("female_allo_17A7","female_allo_17AB","female_allo_17C2","female_allo_17E2","female_camF1","female_Z23698","female_Z23702","female_Z23721","female_Z23726","female_Z23732","female_Z23733","male_allo_17A1","male_allo_1873","male_allo_1876","male_allo_18A1","male_allo_18AF","male_allo_190F","male_camM5","male_Z23701","male_Z23738","male_Z23739")



my_df_wide <- data.frame(contig=character(),
                    start=integer(),
                    stop=integer())


# loop through each sample and make a long format df
#for (sample in sample_vector){
#   sample <- "female_allo_17AB"
#    a <- read.table(paste(eval(sample),"_sorted_rg.bam_RPKMnorm_bamCoverage_10kb.bw", sep=""))
#    a$sample <- eval(sample)
#    colnames(a) <- c("contig","start","stop","depth","sample")
#    my_df <- rbind(my_df,a)
#}  

#head(my_df)
#dim(my_df)

# loop through each sample and correct merged rows
for (sample in sample_vector){
  # sample <- "female_allo_17A7"
  #sample <- "male_allo_18AF"
  a <- read.table(paste(eval(sample),"_sorted_rg.bam_RPKMnorm_bamCoverage_10kb.bw", sep=""))
  #a$sample <- eval(sample)
  colnames(a) <- c("contig","start","stop",paste(eval(sample),"_depth", sep=""))
  # here need to find places where consecutive bins have the same number of reads overlapping, because they will be merged.
  # for example: male_allo_18AF_sorted_rg.bam_RPKMnorm_bamCoverage_10kb.bw:tig00000003	20000	40000	0.0394696
  # Create an empty tibble to store the new data
  new_df <- tibble()
  expected_diff <- 10000
# Loop through each row of the original dataframe
  for (i in 1:nrow(a)) {
    # Add the current row to the new dataframe
    # Check the difference between columns for the current row
    current_diff <- a$stop[i] - a$start[i]
    print(paste(eval(sample)," ",eval(i)," ",eval(current_diff),sep=""))
      # If the difference is not 10000, add a new row below
      if (current_diff != expected_diff) { # there is a row with merged bins
        for (j in 1:ceiling(current_diff/expected_diff)) { #for the number of merged bins
          extra_row <- a[i,] # make a row
          extra_row$start <- a$start[i] + (expected_diff*j) - expected_diff # adjust the start
          extra_row$stop <- a$start[i] + (expected_diff*j)  # adjust the stop
          new_df <- bind_rows(new_df, extra_row) #add this to a new df
        }  
      }else{ # there is not a row with merged bins
          extra_row <- a[i,]
          new_df <- bind_rows(new_df, extra_row)
      }
  }
 # write to a new file
 write.table(new_df, file = paste(eval(sample),"_sorted_rg.bam_RPKMnorm_bamCoverage_10kb_fixedrowz.bw", sep=""), sep = "\t", row.names = FALSE) 
##  my_df_wide <- merge(my_df_wide, new_df, by=c("contig","start","stop"), all=T) # merge the corrected df
  #my_df_wide <- full_join(my_df_wide, new_df, by = c("contig","start","stop")) # merge the corrected df
#  my_df_wide <- merge(my_df_wide, a, by=c("contig","start","stop"), all=T) # merge the corrected df
  # check tig00008642  
}
  
  
  
dim(my_df_wide)

# Now read in the fixed files
# loop through each sample and make a wide format df
for (sample in sample_vector){
  # sample <- "female_allo_17A7"
  #sample <- "male_allo_18AF"
  a <- read.table(paste(eval(sample),"_sorted_rg.bam_RPKMnorm_bamCoverage_10kb_fixedrowz.bw.txt", sep=""), header=T)
  #a$sample <- eval(sample)
  #colnames(a) <- c("contig","start","stop",paste(eval(sample),"_depth", sep=""))
  my_df_wide <- merge(my_df_wide, a, by=c("contig","start","stop"), all=T) # merge the corrected df
}



head(my_df_wide)
dim(my_df_wide)

#library(dplyr)

# use medians instead of means to reduce the effect of outliers
# rowMedians(x, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(x),
#           ..., useNames = TRUE)

# all ----
# calculate the F-M difference in mean depth 
my_df_wide$mean_diff_F_minus_M <- rowMeans(my_df_wide[ , c(4:18)], na.rm=TRUE) -
  rowMeans(my_df_wide[ , c(19:32)], na.rm=TRUE)
my_df_wide$mean_femalez <- rowMeans(my_df_wide[ , c(4:18)], na.rm=TRUE)
my_df_wide$mean_malez <- rowMeans(my_df_wide[ , c(19:32)], na.rm=TRUE) 

# Or only mainland ----
my_df_wide$mean_diff_F_minus_M_onlymainland <- rowMeans(my_df_wide[ , c(4:7,12:18)], na.rm=TRUE) -
  rowMeans(my_df_wide[ , c(19:24,29:32)], na.rm=TRUE)
my_df_wide$mean_femalez_onlymainland <- rowMeans(my_df_wide[ , c(4:7,12:18)], na.rm=TRUE)
my_df_wide$mean_malez_onlymainland <- rowMeans(my_df_wide[ , c(19:24,29:32)], na.rm=TRUE) 

# Or only Bioko ----
my_df_wide$mean_diff_F_minus_M_onlybioko <- rowMeans(my_df_wide[ , c(8:11)], na.rm=TRUE) -
  rowMeans(my_df_wide[ , c(25:28)], na.rm=TRUE)
my_df_wide$mean_femalez_onlybioko <- rowMeans(my_df_wide[ , c(8:11)], na.rm=TRUE)
my_df_wide$mean_malez_onlybioko <- rowMeans(my_df_wide[ , c(25:28)], na.rm=TRUE) 

# also do medians
library(matrixStats)
# all ----
# calculate the F-M difference in mean depth 
my_df_wide$median_diff_F_minus_M <- rowMedians(as.matrix(my_df_wide[ , c(4:18)]), na.rm=TRUE) -
  rowMedians(as.matrix(my_df_wide[ , c(19:32)]), na.rm=TRUE)
my_df_wide$mean_femalez <- rowMedians(as.matrix(my_df_wide[ , c(4:18)]), na.rm=TRUE)
my_df_wide$mean_malez <- rowMedians(as.matrix(my_df_wide[ , c(19:32)]), na.rm=TRUE) 

# Or only mainland ----
my_df_wide$median_diff_F_minus_M_onlymainland <- rowMedians(as.matrix(my_df_wide[ , c(4:7,12:18)]), na.rm=TRUE) -
  rowMedians(as.matrix(my_df_wide[ , c(19:24,29:32)]), na.rm=TRUE)
my_df_wide$median_femalez_onlymainland <- rowMedians(as.matrix(my_df_wide[ , c(4:7,12:18)]), na.rm=TRUE)
my_df_wide$median_malez_onlymainland <- rowMedians(as.matrix(my_df_wide[ , c(19:24,29:32)]), na.rm=TRUE) 

# Or only Bioko ----
my_df_wide$median_diff_F_minus_M_onlybioko <- rowMedians(as.matrix(my_df_wide[ , c(8:11)]), na.rm=TRUE) -
  rowMedians(as.matrix(my_df_wide[ , c(25:28)]), na.rm=TRUE)
my_df_wide$median_femalez_onlybioko <- rowMedians(as.matrix(my_df_wide[ , c(8:11)]), na.rm=TRUE)
my_df_wide$median_malez_onlybioko <- rowMedians(as.matrix(my_df_wide[ , c(25:28)]), na.rm=TRUE) 

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
                         xlab="Chromosome", ylab="Mean Depth Diff (F-M)",
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
    #maxy<-ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0)))+.1;
    #miny<-floor(min(y, ifelse(!is.na(sig.level), sig.level, 0)))-.1;
    maxy<-max(my_df_wide$mean_diff_F_minus_M, na.rm=T)+1;
    miny<-min(my_df_wide$mean_diff_F_minus_M, na.rm=T)-1;
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
                          xlab="Mb", ylab="Mean depth diff (F-M)",
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
                                  xlab="Mb", ylab="Mean Depth Diff (F-M)",
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

# get rid of tig00007285 because this is mtDNA ----
my_df_wide <- my_df_wide[my_df_wide$contig != "tig00007285",]
all_of_em_noNAs <- my_df_wide[complete.cases(my_df_wide), ]

dim(all_of_em_noNAs)
dim(my_df_wide)
# change contig names
#all_of_em_noNAs <- as.data.frame(sapply(all_of_em_noNAs,function(x) {x <- gsub("tig0000","c",x)}))
#all_of_em_noNAs <- as.data.frame(sapply(all_of_em_noNAs,function(x) {x <- gsub("tig000","c",x)}))
# get variables in proper format
#all_of_em_noNAs$start <- as.numeric(all_of_em_noNAs$start)
#all_of_em_noNAs$mean_diff_M_minus_F <- as.numeric(all_of_em_noNAs$mean_diff_M_minus_F)
#all_of_em_noNAs$mean_femalez <- as.numeric(all_of_em_noNAs$mean_femalez)
#all_of_em_noNAs$mean_malez <- as.numeric(all_of_em_noNAs$mean_malez)


# plot all mean diff ----
F_minus_M_allo_meandiff_depth_plot <- manhattan.plot(factor(my_df_wide$contig), my_df_wide$start, 
                                            my_df_wide$mean_diff_F_minus_M)#,

# plot all mean diff ----
F_minus_M_allo_mediandiff_depth_plot <- manhattan.plot(factor(my_df_wide$contig), my_df_wide$start, 
                                            my_df_wide$median_diff_F_minus_M)#,



jpeg("./F_minus_M_allo_meandiff_depth_plot_10kb.jpg",w=300, h=3.0, units ="in", bg="transparent", res = 200)
  F_minus_M_allo_meandiff_depth_plot
dev.off()

jpeg("./F_minus_M_allo_mediandiff_depth_plot_10kb.jpg",w=300, h=3.0, units ="in", bg="transparent", res = 200)
  F_minus_M_allo_mediandiff_depth_plot
dev.off()


View(my_df_wide)

# plot onlymainland ----
F_minus_M_allo_depth_plot_onlymainland <- manhattan.plot(factor(my_df_wide$contig), my_df_wide$start, 
                                            my_df_wide$mean_diff_F_minus_M_onlymainland)#,
jpeg("./F_minus_M_allo_depth_plot_withoutliers_10kb_onlymainland.jpg",w=300, h=3.0, units ="in", bg="transparent", res = 200)
  F_minus_M_allo_depth_plot_onlymainland
dev.off()

# plot onlybioko ----
F_minus_M_allo_depth_plot_onlybioko <- manhattan.plot(factor(my_df_wide$contig), my_df_wide$start, 
                                                         my_df_wide$mean_diff_F_minus_M_onlybioko)#,
jpeg("./F_minus_M_allo_depth_plot_withoutliers_10kb_onlybioko.jpg",w=300, h=3.0, units ="in", bg="transparent", res = 200)
  F_minus_M_allo_depth_plot_onlybioko
dev.off()


# Scatter of mainland vs bioko
library(ggplot2)

# Create a basic scatter plot
y <- ggplot(data = my_df_wide, aes(x = median_diff_F_minus_M_onlymainland, y = median_diff_F_minus_M_onlybioko)) +
  geom_point() +
  theme_classic(base_size = 16) 
jpeg("./mainland_vs_bioko.jpg",w=5, h=5.0, units ="in", bg="transparent", res = 200)
  y
dev.off()


View(my_df_wide)
min(my_df_wide$mean_diff_F_minus_M, na.rm=T)
max(my_df_wide$mean_diff_F_minus_M, na.rm=T)

# Notez
# tig00009267_80000_90000 - no match in XL - high F-M in mainland but not bioko
# tig00011059_7300000_7310000 is 18S rDNA - high F-M in bioko but not mainland

# tig00008522:3580000-3600000 -zc3h6 - low F-M in all # no known sex determining function


my_df_wide_onlymalecoverage <- my_df_wide[my_df_wide$mean_femalez == 0,]
# get rid of rows with no data
my_df_wide_onlymalecoverage <- my_df_wide_onlymalecoverage[rowSums(is.na(my_df_wide_onlymalecoverage)) != ncol(my_df_wide_onlymalecoverage), ]
dim(my_df_wide_onlymalecoverage)
View(my_df_wide_onlymalecoverage)


# plot ----
M_minus_F_allo_depth_plot <- manhattan.plot(factor(all_of_em_noNAs$contig), all_of_em_noNAs$start, 
                                            all_of_em_noNAs$mean_diff_M_minus_F)#,
                                          # ylim= c(0,12),xlab = "",
                                         # par.settings = theme.novpadding)#,
                                          #panel=function(x, y, ...){
                                            # panel.rect(xleft=695, ybottom=0,
                                            #           xright=733, ytop=10, alpha=0.3, col="light blue")
                                          #panel.text(275,9,labels=expression(italic("X. allofraseri")),fontsize=14)})

jpeg("./M_minus_F_allo_depth_plot_withoutliers_10kb.jpg",w=300, h=3.0, units ="in", bg="transparent", res = 200)
  M_minus_F_allo_depth_plot
dev.off()

# Weirdonez
SexChr <- my_df_wide[my_df_wide$contig == "tig00011059",] # this includes 18S rDNA
SexChr <- my_df_wide[my_df_wide$contig == "tig00008640",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00009533",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00006777",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00011479",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00008522",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00000697",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00010978",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00005103",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00007523",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00007559",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00008322",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00008988",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00000216",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00004247",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00004787",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00004937",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00005278",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00007882",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00009681",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00010342",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00010588",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00010385",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00010638",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00010774",]

SexChr <- my_df_wide[my_df_wide$contig == "tig00011059",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00011782",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00011811",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00009379",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00008278",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00009515",]

SexChr <- my_df_wide[my_df_wide$contig == "tig00005092",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00003683",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00008584",]

SexChr <- my_df_wide[my_df_wide$contig == "tig00005199",]
SexChr <- my_df_wide[my_df_wide$contig == "tig00009996",]

SexChr <- my_df_wide[my_df_wide$contig == "tig00008985",]


# narrow down region
# SexChr <- SexChr[((SexChr$start > 140000)&(SexChr$start < 280000)),]
# SexChr <- SexChr[(SexChr$start < 200000),]
Contig_focus_plot <- chromosome.plot.no.thin(factor(SexChr$contig),SexChr$start,SexChr$median_diff_F_minus_M)
Contig_focus_plot <- chromosome.plot.no.thin(factor(SexChr$contig),SexChr$start,SexChr$median_diff_F_minus_M_onlymainland);Contig_focus_plot
Contig_focus_plot <- chromosome.plot.no.thin(factor(SexChr$contig),SexChr$start,SexChr$median_diff_F_minus_M_onlybioko);Contig_focus_plot
                                              

jpeg("./tig00009996_mediandif_10kbfocus.jpg",w=7, h=3.0, units ="in", bg="transparent", res = 200)
  Contig_focus_plot
dev.off()



View(SexChr)

```
