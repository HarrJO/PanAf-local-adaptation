# baypass_tools.R
## General tools for analysing the output from the BayPass cor or AUX model
# Library
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(shadowtext)
library(viridis)
library(ggrepel)
library(ggspatial)

# Plot
manhattan_plot=function(df, title="Manhattan Plot", stat="`BF(dB).median`", chr.lengths=NULL, gtf=NULL, focal.snps.list=NULL, focal.snps.cols=c("purple", "red", "darkorange"),
                        cov_col=NULL, xlim=NULL, gtf_features=c("gene"), gene_names=FALSE){
  if(stat %in% c("`BF(dB)`", "`BF(dB).median`", "BF.med")){
    stat_lab="BF(dB)"
  }else if(stat %in% c("XtXst", "XtXst.med")){
    stat_lab="XtX*"
  }else if(stat=='log_p'){
    stat_lab=expression(-log[10]~empirical~p~value)
  }else if(stat=='log_fpr'){
    stat_lab=expression(-log[10]~FPR)
  }else{stat_lab=stat}
  df=df[order(df$chr),]
  df$chr=as.factor(df$chr)
  chrs=unique(df$chr)
  # Prepare cumulative position column for plotting
  ## Hg19 chr lengths by default
  if(is.null(chr.lengths)){
    chr.lengths=data.frame(chr=1:22, length=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,
                                              133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566))
  }
  ### Only include selected chrs
  chr.lengths=chr.lengths[chr.lengths$chr %in% chrs,]
  ## Get cumulative sums of chr lengths
  cumsum=cumsum(chr.lengths$length)
  ## Get vector of chr starts
  if(length(chrs)>1){
    cum.chr.start=data.frame(chr=chrs, cum.chr.start=c(0, cumsum[1:(length(cumsum)-1)]))
  }else{
    cum.chr.start=data.frame(chr=chrs, cum.chr.start=0)
  }
  ## Add cumulative position to df data
  df=merge(df, cum.chr.start, by='chr')
  df$cumpos=as.numeric(df$pos + df$cum.chr.start)
  # Prepare table for plotting the positions of genes if a gtf file is provided
  if(!is.null(gtf)){
    ## Keep only selected chrs
    gtf=gtf[gtf$V1 %in% chrs,]
    ## Get the maximum and minimum positions of the gene (as the gtf may contain CDS and exon information instead of whole gene coordinates)
    ### V10 is the gene name column, V1 is chr, V4 is start and V5 is end
    #gene.bed=gtf %>% group_by(V10) %>% summarize(chr=unique(V1), start=min(V4), end=max(V5))
    #colnames(gene.bed)[1]='gene'
    ## Add cumulative position of chr
    #gene.bed=merge(gene.bed, cum.chr.start, by='chr')
    ## Make cumulative start and ends for gene
    #gene.bed$cumstart=gene.bed$start + gene.bed$cum.chr.start
    #gene.bed$cumend=gene.bed$end + gene.bed$cum.chr.start
    #gene.bed$cummid=(gene.bed$cumstart+gene.bed$cumend)/2
    
    gtf=gtf[gtf$V3 %in% gtf_features,]
    gtf$gene=paste(gtf$V3, gtf$V10, 1:nrow(gtf),sep="_")
    #gtf$gene_name=paste(gtf$V3, gtf$V10, sep="-")
    gtf$gene_name=gtf$V10
    
    #starts=aggregate(V4~V1+gene, data = gtf, FUN = min) # Old version of R
    starts=aggregate(x = gtf$V4, by = list(gtf$V1, gtf$gene), FUN = min) # For new version of R
    colnames(starts)=c('chr', 'gene', 'start')
    #ends=aggregate(V5~V1+gene, data = gtf, FUN = max) # Old version of R
    ends=aggregate(x = gtf$V5, by = list(gtf$V1, gtf$gene), FUN = max) # For new version of R
    colnames(ends)=c('chr', 'gene', 'end')
    gene.bed=merge(starts, ends)
    
    ## Add cumulative position of chr
    gene.bed=merge(gene.bed, cum.chr.start, by='chr')
    ## Add gene names
    gene.bed=merge(gene.bed, unique(gtf[, c('gene', 'gene_name')]), by='gene')
    ## Make cumulative start and ends for gene
    gene.bed$cumstart=gene.bed$start + gene.bed$cum.chr.start
    gene.bed$cumend=gene.bed$end + gene.bed$cum.chr.start
    gene.bed$cummid=(gene.bed$cumstart+gene.bed$cumend)/2
  }
  # Manhatten Plot
  ## Get thresholds
  stat.thresh=c()
  ## Prepare x axis
  ### Midpoints
  chr.mid=data.frame(chrs=chrs, center=(cumsum+cum.chr.start$cum.chr.start)/2)
  ### Min and Max
  chr.max=cumsum
  chr.min=cum.chr.start$cum.chr.start
  ## Plotting
  ### Set plotting colours
  chr.panel.col=rep(c("grey", "white"), 20)
  #### Crop to the number of chrs
  chr.panel.col=chr.panel.col[1:length(chrs)]
  #### lims
  min.y=min(df[[stat]])
  max.y=max(df[[stat]])
  range.y=max.y-min.y
  min.y=min.y-0.1*range.y
  max.y=max.y+0.05*range.y
  ### Plot
  if(length(chrs)==1){
    xlab=paste0('Chromosome ', chrs)
  }else{
    xlab='Chromosome'
    }
  man.plot=ggplot(NULL) +
    #### General
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour = 'black', face="bold", size=7)) +
    labs(title=paste0(title), x=xlab, y=stat_lab) +
    #### Axes
    scale_x_continuous(label = chr.mid$chr, breaks= chr.mid$center ) +
    scale_y_continuous(limits = c(min.y, max.y), expand = c(0, 0) ) +
    geom_rect(data=NULL, aes(xmin=chr.min, xmax=chr.max, ymin=-Inf, max=Inf), fill=chr.panel.col, alpha=0.3)
  #### Get max and min y
  #min.y=layer_scales(man.plot)$y$range$range[2]
  #max.y=layer_scales(man.plot)$y$range$range[0]
  #### Add manual colour scale?
  if(!is.null(cov_col)){
    man.plot=man.plot+scale_discrete_manual(name='Covariate', values = cov_col, aesthetics = c("colour", "fill"))
  }
  #### Add gene positions?
  if(!is.null(gtf)){
    man.plot=man.plot+geom_rect(data=gene.bed, aes(xmin=cumstart, xmax=cumend, ymin=-Inf, max=Inf, fill=gene), alpha=0.2)+ 
      scale_fill_discrete(guide="none")
    #### Add gene names?
    if(gene_names){
      man.plot=man.plot+
        geom_shadowtext(data=gene.bed, aes(x=cummid, y=min.y+0.05*(max.y-min.y), label=gene_name, colour=gene), alpha=1, bg.colour='black') +
        scale_colour_discrete(guide = "none")
    }
  }
  # Points
  if(nrow(df)>1000){
    point_alpha=0.1
    point_size=1
    focal_point_size=1
  }else{
    point_alpha=1
    point_size=1
    focal_point_size=2
  }
  man.plot=man.plot+
    geom_point(data=df, aes_string(x='cumpos', y=stat, col="COVARIABLE_name"), alpha=point_alpha, size=1, shape=18, col="turquoise4")
    
  ### Highlight SNPs?
  if(!is.null(focal.snps.list)){
    for(i in 1:length(focal.snps.list)){
      focal.snps=focal.snps.list[[i]]
      focal.snps$chr=as.factor(focal.snps$chr)
      ### Get the results for focal SNPs (doing a left merge (i.e. all.x=TRUE) we only keep results for our focal SNPs)
      focal.snps=merge(focal.snps[,c('chr','pos')], df, by=c('chr', 'pos'), all.x=TRUE)
      ### Add to plot
      man.plot=man.plot+geom_point(data=focal.snps, aes_string(x='cumpos', y=stat), colour=focal.snps.cols[i], alpha=1, size=focal_point_size, shape=18)
    }
  }
  #### xlim?
  if(!is.null(xlim)){
    if(length(chrs)==1){
      man.plot=man.plot+xlim(xlim)
    }else{
      stop("Cannot zoom into this region, multiple chromosomes present")
    }
  }
  # If there are too many covariables, plot the legend separately
  n_cat=length(unique(df$COVARIABLE_name))
  if(n_cat>10){
    ## Extract legend
    legend=get_legend(man.plot)
    ## Remove legend
    man.plot=man.plot+theme(legend.position = "none")
  }
  # If there is only one covariable, dont bother to plot the legend 
  if(n_cat==1){
    ## Remove legend
    man.plot=man.plot+theme(legend.position = "none")
  }
  print(man.plot)
  ## Plot legend
  if(n_cat>10){
    grid.newpage()
    grid.draw(legend)
  }
}

manhattan_plot.v2=function(df, title="Manhattan Plot", stat="`BF(dB).median`", chr.lengths=NULL, gtf=NULL, focal.snps.list=NULL, focal.snps.cols=c("purple", "red", "darkorange"),
                        cov_col=NULL, xlim=NULL, gtf_features=c("gene"), gene_names=FALSE, trans_dir=FALSE, line=FALSE){
  if(stat %in% c("`BF(dB)`", "`BF(dB).median`", "BF.med")){
    stat_lab="BF(dB)"
  }else if(stat %in% c("XtXst", "XtXst.med")){
    stat_lab="XtX*"
  }else if(stat=='log_p'){
    stat_lab=expression(-log[10]~empirical~p~value)
  }else if(stat=='log_fpr'){
    stat_lab=expression(-log[10]~FPR)
  }else{stat_lab=stat}
  df=df[order(df$chr),]
  df$chr=as.factor(df$chr)
  chrs=unique(df$chr)
  # Prepare cumulative position column for plotting
  ## Hg19 chr lengths by default
  if(is.null(chr.lengths)){
    chr.lengths=data.frame(chr=1:22, length=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,
                                              133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566))
  }
  ### Only include selected chrs
  chr.lengths=chr.lengths[chr.lengths$chr %in% chrs,]
  ## Get cumulative sums of chr lengths
  cumsum=cumsum(chr.lengths$length)
  ## Get vector of chr starts
  if(length(chrs)>1){
    cum.chr.start=data.frame(chr=chrs, cum.chr.start=c(0, cumsum[1:(length(cumsum)-1)]))
  }else{
    cum.chr.start=data.frame(chr=chrs, cum.chr.start=0)
  }
  ## Add cumulative position to df data
  df=merge(df, cum.chr.start, by='chr')
  df$cumpos=as.numeric(df$pos + df$cum.chr.start)
  # Prepare table for plotting the positions of genes if a gtf file is provided
  if(!is.null(gtf)){
    ## Keep only selected chrs
    gtf=gtf[gtf$V1 %in% chrs,]
    ## Get the maximum and minimum positions of the gene (as the gtf may contain CDS and exon information instead of whole gene coordinates)
    ### V10 is the gene name column, V1 is chr, V4 is start and V5 is end
    #bed=gtf %>% group_by(V10) %>% summarize(chr=unique(V1), start=min(V4), end=max(V5))
    #colnames(bed)[1]='gene'
    ## Add cumulative position of chr
    #bed=merge(bed, cum.chr.start, by='chr')
    ## Make cumulative start and ends for gene
    #bed$cumstart=bed$start + bed$cum.chr.start
    #bed$cumend=bed$end + bed$cum.chr.start
    #bed$cummid=(bed$cumstart+bed$cumend)/2
    
    gtf=gtf[gtf$V3 %in% c("gene", gtf_features),]
    gtf$gene=paste(gtf$V3, gtf$V10, 1:nrow(gtf),sep="_")
    #gtf$gene_name=paste(gtf$V3, gtf$V10, sep="-")
    gtf$gene_name=gtf$V10
    
    #starts=aggregate(V4~V1+V7+gene, data = gtf, FUN = min)#  Old version of R
    starts=aggregate(x = gtf$V4, by = list(gtf$V1, gtf$V7, gtf$gene), FUN = min) # For new version of R
    colnames(starts)=c('chr', 'strand', 'gene', 'start')
    #ends=aggregate(V5~V1+V7+gene, data = gtf, FUN = max)# Old version of R
    ends=aggregate(x = gtf$V5, by = list(gtf$V1, gtf$V7, gtf$gene), FUN = max) # For new version of R
    colnames(ends)=c('chr', 'strand','gene', 'end')
    bed=merge(starts, ends)
    
    ## Add cumulative position of chr
    bed=merge(bed, cum.chr.start, by='chr')
    ## Add gene names
    bed=merge(bed, unique(gtf[, c('gene', 'gene_name')]), by='gene')
    ## Make cumulative start and ends for gene
    bed$cumstart=bed$start + bed$cum.chr.start
    bed$cumend=bed$end + bed$cum.chr.start
    bed$cummid=(bed$cumstart+bed$cumend)/2
    
    # Gene gene mid
    #gene.mid=aggregate(formula = cummid ~ gene_name, data = bed, FUN=mean) # Old version of R
    gene.mid=aggregate(x = bed$cummid, by = list(bed$gene_name), FUN = mean) # New version of R
    colnames(gene.mid)=c('gene_name', 'cummid.mean')
    bed=merge(bed, gene.mid, by='gene_name')
    
    # Seperate gene df and selected feature df
    ## We need seperate ones because even if we are plotting exons we want the labels to be gene names
    gene.bed=bed[grepl(paste0("^gene_"), bed$gene), ]
    select.bed=bed[grepl(paste(paste0("^",gtf_features,"_"), collapse="|"), bed$gene), ]
    
  }
  # Manhatten Plot
  ## Get thresholds
  stat.thresh=c()
  ## Prepare x axis
  ### Midpoints
  chr.mid=data.frame(chrs=chrs, center=(cumsum+cum.chr.start$cum.chr.start)/2)
  ### Min and Max
  chr.max=cumsum
  chr.min=cum.chr.start$cum.chr.start
  ## Plotting
  ### Set plotting colours
  chr.panel.col=rep(c("grey", "white"), 20)
  #### Crop to the number of chrs
  chr.panel.col=chr.panel.col[1:length(chrs)]
  #### lims
  min.y=min(df[[stat]])
  max.y=max(df[[stat]])
  range.y=max.y-min.y
  min.y=min.y-0.4*range.y
  max.y=max.y+0.05*range.y
  ### Plot
  if(length(chrs)==1){
    xlab=paste0('Chromosome ', chrs)
  }else{
    xlab='Chromosome'
  }
  man.plot=ggplot(NULL) +
    #### General
    theme_bw() +
    theme(
      panel.border = element_blank(),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour = 'black', face="bold", size=7)) +
    labs(title=paste0(title), x=xlab, y=stat_lab) +
    #### Axes
    scale_x_continuous(label = chr.mid$chr, breaks= chr.mid$center ) +
    scale_y_continuous(limits = c(min.y, max.y), expand = c(0, 0) ) +
    geom_rect(data=NULL, aes(xmin=chr.min, xmax=chr.max, ymin=-Inf, max=Inf), fill=chr.panel.col, alpha=0.3)
  #### Get max and min y
  #min.y=layer_scales(man.plot)$y$range$range[2]
  #max.y=layer_scales(man.plot)$y$range$range[0]
  #### Add manual colour scale?
  if(!is.null(cov_col)){
    man.plot=man.plot+scale_discrete_manual(name='Covariate', values = cov_col, aesthetics = c("colour", "fill"))
  }
  #### Add gene positions?
  if(!is.null(gtf)){
    man.plot=man.plot+
      geom_hline(yintercept = min.y+(0.275*range.y)) +
      geom_rect(data=select.bed, aes(xmin=cumstart, xmax=cumend, ymin=min.y+(0.25*range.y), ymax=min.y+(0.3*range.y), fill=gene_name, col=gene_name), alpha=1)+ 
      scale_fill_discrete(guide="none")
    #### Add gene names?
    if(gene_names){
      min_break=min(df[[stat]], rm.na=T)
      max_break=max(df[[stat]], rm.na=T)
      breaks=round(seq(min_break, max_break, (max_break-min_break)/3), 2)
      #breaks=c(min_break, mean(c(min_break, max_break)), max_break)
      man.plot=man.plot+
        geom_shadowtext(data=gene.bed, aes(x=cummid.mean, y=min.y+(0.05*range.y), label=gene_name, colour=gene_name), alpha=1, bg.colour='black', size=3) +
        guides(colour = FALSE, fill=FALSE) +
        scale_colour_discrete(guide = "none") +
        scale_y_continuous(breaks=breaks)
        #theme(axis.text.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
    }
    #### Add transcription direction?
    if(trans_dir){
      gene.bed.pos=gene.bed[gene.bed$strand=="+",]
      gene.bed.neg=gene.bed[gene.bed$strand=="-",]
      int=0.05/(nrow(gene.bed)-1)
      # If there is more than one transcript
      if(int>0){
        y=min.y+seq(0.15, 0.2, int)*range.y
      }else{
        y=0.175
      }
      names(y)=gene.bed$gene_name
      if(nrow(gene.bed.pos)>0){
        man.plot=man.plot+
          geom_segment(data=gene.bed.pos, aes(x = cumstart, xend = cumend, y = y[gene.bed.pos$gene_name], yend=y[gene.bed.pos$gene_name], colour=gene_name),
                       arrow = arrow(length = unit(0.1, "cm")), size = 0.5)
      }
      if(nrow(gene.bed.neg)>0){
        man.plot=man.plot+
          geom_segment(data=gene.bed.neg, aes(x = cumend, xend = cumstart, y = y[gene.bed.neg$gene_name], yend=y[gene.bed.neg$gene_name], colour=gene_name),
                       arrow = arrow(length = unit(0.1, "cm")), size = 0.5)
      }
    }
    colours=rainbow(length(unique(select.bed[['gene_name']])))
    names(colours)=unique(select.bed[['gene_name']])
    man.plot=man.plot+
      #scale_discrete_manual(values=colorRampPalette(brewer.pal(8, "Set2"))(length(unique(select.bed[['gene_name']]))), aesthetics = c("colour", "fill")) 
      scale_discrete_manual(values=colours, aesthetics = c("colour", "fill")) 
  }
  # Points
  if(nrow(df)>1000){
    point_alpha=0.1
    point_size=1
    focal_point_size=1
  }else{
    point_alpha=1
    point_size=1
    focal_point_size=2
  }
  man.plot=man.plot+
    geom_point(data=df, aes_string(x='cumpos', y=stat, col="COVARIABLE_name"), alpha=point_alpha, size=1, shape=18, col="turquoise4") #+
    #scale_discrete_manual(values= wes_palette("GrandBudapest1", n = length(unique(select.bed[['gene_name']]))), aesthetics = c("colour", "fill")) 
    #scale_discrete_manual(values= rainbow(n = length(unique(select.bed[['gene_name']]))), aesthetics = c("colour", "fill")) 
    #scale_discrete_manual(values=brewer.pal(n = length(unique(select.bed[['gene_name']])), name = "Set1"), aesthetics = c("colour", "fill")) 
  ### Highlight SNPs?
  if(!is.null(focal.snps.list)){
    for(i in 1:length(focal.snps.list)){
      focal.snps=focal.snps.list[[i]]
      focal.snps$chr=as.factor(focal.snps$chr)
      ### Get the results for focal SNPs (doing a left merge (i.e. all.x=TRUE) we only keep results for our focal SNPs)
      focal.snps=merge(focal.snps[,c('chr','pos')], df, by=c('chr', 'pos'), all.x=TRUE)
      ### Add to plot
      man.plot=man.plot+geom_point(data=focal.snps, aes_string(x='cumpos', y=stat), colour=focal.snps.cols[i], alpha=1, size=focal_point_size, shape=18)
    }
  }
  # line
  if(line){
    man.plot=man.plot+
      geom_smooth(data=df, aes_string(x='cumpos', y=stat, col="COVARIABLE_name"), method = "loess", col="turquoise4")
  }
  #### xlim?
  if(!is.null(xlim)){
    if(length(chrs)==1){
      man.plot=man.plot+xlim(xlim)
    }else{
      stop("Cannot zoom into this region, multiple chromosomes present")
    }
  }
  # If there are too many covariables, plot the legend separately
  n_cat=length(unique(df$COVARIABLE_name))
  if(n_cat>10){
    ## Extract legend
    legend=get_legend(man.plot)
    ## Remove legend
    man.plot=man.plot+theme(legend.position = "none")
  }
  # If there is only one covariable, dont bother to plot the legend 
  if(n_cat==1){
    ## Remove legend
    man.plot=man.plot+theme(legend.position = "none")
  }
  print(man.plot)
  ## Plot legend
  if(n_cat>10){
    grid.newpage()
    grid.draw(legend)
  }
}

manhattan_plot_per_gene=function(subsps, genes, gtf_features="exon"){
  for(subsp in subsps){
    for(gene in genes){
      gtf.tmp=gtf[gtf$V10 %in% gene,]
      betai_fpr.tmp=betai_fpr[[subsp]][betai_fpr[[subsp]]$gene %in% gene,]
      manhattan_plot(betai_fpr.tmp, stat="log_fpr", cov_col="black", title=paste(subsp, "; ", gene, "; ", gtf_features), focal.snps.list=top_2[subsp],
                     gtf=gtf.tmp, xlim=c(min(gtf.tmp$V4), max(gtf.tmp$V5)), gtf_features=gtf_features, focal.snps.cols=c("red"))
      cat(gene)
      print(head(betai_fpr.tmp[order(betai_fpr.tmp$fpr),]))
    }
  }
}

snp_gene_positions=function(df, title="", chr.lengths=NULL, gtf=NULL, gtf_features=c("gene"), xlim=NULL){
  df=df[order(df$chr),]
  df$chr=as.factor(df$chr)
  chrs=unique(df$chr)
  # Prepare cumulative position column for plotting
  ## Hg19 chr lengths by default
  if(is.null(chr.lengths)){
    chr.lengths=data.frame(chr=1:22, length=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,
                                              133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566))
  }
  ### Only include selected chrs
  chr.lengths=chr.lengths[chr.lengths$chr %in% chrs,]
  ## Get cumulative sums of chr lengths
  cumsum=cumsum(chr.lengths$length)
  ## Get vector of chr starts
  if(length(chrs)>1){
    cum.chr.start=data.frame(chr=chrs, cum.chr.start=c(0, cumsum[1:(length(cumsum)-1)]))
  }else{
    cum.chr.start=data.frame(chr=chrs, cum.chr.start=0)
  }
  ## Add cumulative position to df data
  df=merge(df, cum.chr.start, by='chr')
  df$cumpos=as.numeric(df$pos + df$cum.chr.start)
  # Prepare table for plotting the positions of genes if a gtf file is provided
  if(!is.null(gtf)){
    ## Keep only selected chrs
    gtf=gtf[gtf$V1 %in% chrs,]
    ## Get the maximum and minimum positions of the gene (as the gtf may contain CDS and exon information instead of whole gene coordinates)
    ### V10 is the gene name column, V1 is chr, V4 is start and V5 is end
    #gene.bed=gtf %>% group_by(V10) %>% summarize(chr=unique(V1), start=min(V4), end=max(V5))
    #colnames(gene.bed)[1]='gene'
    ## Add cumulative position of chr
    #gene.bed=merge(gene.bed, cum.chr.start, by='chr')
    ## Make cumulative start and ends for gene
    #gene.bed$cumstart=gene.bed$start + gene.bed$cum.chr.start
    #gene.bed$cumend=gene.bed$end + gene.bed$cum.chr.start
    #gene.bed$cummid=(gene.bed$cumstart+gene.bed$cumend)/2
    
    gtf=gtf[gtf$V3 %in% gtf_features,]
    gtf$gene=paste(gtf$V3, gtf$V10, 1:nrow(gtf),sep="_")
    
    #starts=aggregate(V4~V1+gene, data = gtf, FUN = min) # Old version of R
    starts=aggregate(x = gtf$V4, by = list(gtf$V1, gtf$gene), FUN = min)  # New version of R
    colnames(starts)=c('chr', 'gene', 'start')
    #ends=aggregate(V5~V1+gene, data = gtf, FUN = max) # Old version of R
    ends=aggregate(x = gtf$V5, by = list(gtf$V1, gtf$gene), FUN = max)  # New version of R
    colnames(ends)=c('chr', 'gene', 'end')
    gene.bed=merge(starts, ends)
    
    ## Add cumulative position of chr
    gene.bed=merge(gene.bed, cum.chr.start, by='chr')
    ## Make cumulative start and ends for gene
    gene.bed$cumstart=gene.bed$start + gene.bed$cum.chr.start
    gene.bed$cumend=gene.bed$end + gene.bed$cum.chr.start
    gene.bed$cummid=(gene.bed$cumstart+gene.bed$cumend)/2
  }
  # Plot
  ## Prepare x axis
  ### Midpoints
  chr.mid=data.frame(chrs=chrs, center=(cumsum+cum.chr.start$cum.chr.start)/2)
  ### Min and Max
  chr.max=cumsum
  chr.min=cum.chr.start$cum.chr.start
  ## Plotting
  ### Set plotting colours
  chr.panel.col=rep(c("grey", "white"), 20)
  #### Crop to the number of chrs
  chr.panel.col=chr.panel.col[1:length(chrs)]
  ### Plot
  plot=ggplot(NULL) +
    #### General
    theme_void() +
    theme(panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour = 'black', face="bold", size=7)) +
    labs(title=paste0(title), x="Chromosome") +
    #### Axes
    scale_x_continuous(label = chr.mid$chr, breaks= chr.mid$center ) +
    geom_rect(data=NULL, aes(xmin=chr.min, xmax=chr.max, ymin=-Inf, max=Inf), fill=chr.panel.col, alpha=0.3) +
    geom_point(data=df, aes(x=cumpos, y=runif(nrow(df), min=-0.5, 0.5), col=colour), alpha=1, size=1.5, shape=18) +
    ylim(1,-1)
  #### Add gene positions?
  if(!is.null(gtf)){
    plot=plot+geom_rect(data=gene.bed, aes(xmin=cumstart, xmax=cumend, ymin=-Inf, max=Inf, fill=gene), alpha=0.25)+ 
      scale_fill_discrete(guide="none")
    #### Add gene names if not many
    if(nrow(gene.bed)<20){
      ##### For some reason you need to assign this to an object first then put it into geom_label()
      max.y=layer_scales(plot)$y$range$range[2]
      ##### runif randomly selects a y value in the range  to minimise overlap
      plot=plot+geom_text(data=gene.bed, aes(x=cummid, y=runif(nrow(gene.bed), min=0.4, 0.98), label=gene), alpha=1)
    }
  }
  #### xlim?
  if(!is.null(xlim)){
    if(length(chrs)==1){
      plot=plot+xlim(xlim)
    }else{
      stop("Cannot zoom into this region, multiple chromosomes present")
    }
  }
  print(plot)
}

pairwise_distances=function(cands, bg_snps, subtitle="", reps=30){
  dists=c()
  dists_min=c()
  null_dists=list()
  null_dists_min=list()
  chrs=unique(cands$chr)
  for(chr in chrs[order(chrs)]){
    # Select chromosome
    cands.chr=cands[cands$chr==chr,]
    if(nrow(cands.chr)>1){
      cands.chr=cands.chr[order(cands.chr$pos),]
      bg_snps_chr=bg_snps[bg_snps$chr==chr,]
      bg_snps_chr=bg_snps_chr[order(bg_snps_chr$pos),]
      # Get distance matrix
      dists_chr=dist(cands.chr$pos)
      dists=c(dists, dists_chr)
      # Get minimum distances
      dists_chr_m=as.matrix(dists_chr)
      dists_chr_m[dists_chr_m==0]=NA
      dists_chr_min=apply(dists_chr_m, 1, FUN = min, na.rm = TRUE)
      dists_min=c(dists_min, dists_chr_min)
      # Null SNPs
      for(i in as.character(1:reps)){
        ### Randomly sample the same number of SNPs in the tail from all sites
        null=bg_snps_chr[sample(nrow(bg_snps_chr), nrow(cands.chr)),]
        # Distance of null SNPs
        null_dists_chr=dist(null$pos)
        ## Get minimum distances
        null_dists_chr_m=as.matrix(null_dists_chr)
        null_dists_chr_m[null_dists_chr_m==0]=NA
        null_dists_chr_min=apply(null_dists_chr_m, 1, FUN = min, na.rm = TRUE)
        if(length(null_dists[[i]])>0){
          null_dists[[i]]=c(null_dists[[i]], as.vector(null_dists_chr))
          null_dists_min[[i]]=c(null_dists_min[[i]], as.vector(null_dists_chr_min))
        }else{
          null_dists[[i]]=as.vector(null_dists_chr)
          null_dists_min[[i]]=as.vector(null_dists_chr_min)
        }
      }
    }
  }
  # Plot all pairwise distances
  plot1=ggplot(NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.text=element_text(size=10)) +
    labs(title=paste0("Pairwise Distances Between SNPs"), subtitle=paste0(subtitle), x="Physical Distance (bp)", y = "Density") +
    geom_histogram(aes(x=as.vector(dists), y = ..density.., col="Candidates", fill="Candidates"), alpha=0.1) +
    geom_histogram(aes(x=as.vector(unlist(null_dists)), y = ..density.., col="Total Null", fill="Total Null"), alpha=0.1) +
    scale_discrete_manual(name='Data', values = c('Candidates' = 'red', "Total Null" = 'blue'), 
                          labels=c('Candidates',paste0('Null (total over ',reps,'\nrandom resamples)')),
                          aesthetics = c("colour", "fill"))

  plot2=ggplot(NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.text=element_text(size=10)) +
    labs(title=paste0("Minimum Pairwise Distances Between SNPs"), subtitle=paste0(subtitle), x="Physical Distance (bp)", y = "Density") +
    geom_histogram(aes(x=as.vector(dists_min), y = ..density.., col="Candidates", fill="Candidates"), alpha=0.1) +
    geom_histogram(aes(x=as.vector(unlist(null_dists_min)), y = ..density.., col="Total Null", fill="Total Null"), alpha=0.1) +
    scale_discrete_manual(name='Data', values = c('Candidates' = 'red', "Total Null" = 'blue'), 
                          labels=c('Candidates',paste0('Null (total over ',reps,'\nrandom resamples)')),
                          aesthetics = c("colour", "fill"))
  for(plot in list(plot1, plot2)){print(plot)}
}

sample_map=function(df, title=NULL, col=NULL){
  if(sum(c("Community", "Latitude", "Longitude") %in% colnames(df))<3){
    stop("Input must have at least four columns ; 'Community', 'Latitude', 'Longitude' and then a column for each environmnetal variable to be plotted")
  }
  # Load maps
  world_map=map_data("world")
  chimp.ranges=readOGR("~/OneDrive - University College London/Projects/Ostridge_PanAf/environmental_data/data/maps/chimp_ranges/chimp_ranges.shp")
  ## Downloaded version 5.0.0 from https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-rivers-lake-centerlines/
  rivers=readOGR("~/OneDrive - University College London/Projects/Ostridge_PanAf/environmental_data/data/maps/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp")
  rivers=rivers[grepl("Congo|Niger|Ubangi|Ogooué|Tanganyika|Sanaga|Benue|Uelé|Lualaba|Chambeshi|Chinko|Mbomou|Uele", rivers$name_en),]
  ## Download water bodies Zip (last updated Nov 14, 2018) https://datacatalog.worldbank.org/search/dataset/0040797
  lakes=readOGR("~/OneDrive - University College London/Projects/Ostridge_PanAf/environmental_data/data/maps/africawaterbody/Africa_waterbody.shp")
  # Limits
  ## x
  xmin=min(df$Longitude)
  xmax=max(df$Longitude)
  xmin=xmin-(xmax-xmin)*0.05
  xmax=xmax+(xmax-xmin)*0.05
  ## x
  ymin=min(df$Latitude)
  ymax=max(df$Latitude)
  ymin=ymin-(ymax-ymin)*0.15
  ymax=ymax+(ymax-ymin)*0.15
  # Plot
  map=ggplot(NULL) +
    ## Map
    theme_void() +
    #theme(panel.background = element_rect(fill = 'cadetblue1')) +
    geom_polygon(data=world_map, aes(x = long, y = lat, group = group), fill="navajowhite3", colour = "navajowhite4", alpha=0.25) +
    coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) +
    labs(title=title) +
    theme(plot.title = element_text(hjust = 0.5, size=25)) +
    ## Subspecies ranges
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="troglodytes",], aes(x = long, y = lat, group = group), fill="green3", alpha=0.4) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="schweinfurthii",], aes(x = long, y = lat, group = group), fill="darkorange", alpha=0.4) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="ellioti",], aes(x = long, y = lat, group = group), fill="red", alpha=0.4) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="verus",], aes(x = long, y = lat, group = group), fill="blue", alpha=0.4) +
    ## Water
    geom_polygon(data=lakes, aes(x = long, y = lat, group=group), fill="deepskyblue", alpha=0.5)+
    geom_path(data=rivers, aes(x = long, y = lat, group = group), col="deepskyblue", alpha=0.5) +
    ## Plot points
    geom_point(data=df, aes_string(x='Longitude', y='Latitude'), col='black', fill='black', size=1.5, shape=23) +
    ## Population labels
    geom_text_repel(data=df, aes(x=Longitude, y=Latitude, label=Community), col='black',
                   fontface= "bold",
                   min.segment.length = unit(0, 'lines'),
                   size=3,
                   seed=2,
                   bg.color = "white",
                   bg.r = 0.1,
                   max.overlaps = Inf) 
  print(map)
}

calculate_new_longitude <- function(lon_old, lat_old, distance_km, bearing) {
  # Radius of the Earth (mean radius in kilometers)
  R <- 6371
  # Convert latitude and longitude to radians
  lon_old_rad <- lon_old * pi / 180
  lat_old_rad <- lat_old * pi / 180
  bearing_rad <- bearing * pi / 180  # Convert bearing to radians
  # Calculate new longitude using the haversine formula
  lon_new_rad <- lon_old_rad + atan2(sin(bearing_rad) * sin(distance_km / R) * cos(lat_old_rad),
                                     cos(distance_km / R) - sin(lat_old_rad) * sin(lat_old_rad))
  # Convert back to degrees
  lon_new <- lon_new_rad * 180 / pi
  return(lon_new)
}

var_map=function(df, size_col=NULL, fill_col=NULL, colour=NULL, title=NULL, min_size=2, max_size=6, log_fill=FALSE, seed=6, text_size=2.5, box.padding=0.5, force=1,
                 legend.position=c(0.12, 0.24), ymax_over_ride=NULL, ymin_over_ride=NULL, world_map=NULL, chimp.ranges=NULL, rivers=NULL, lakes=NULL, scale=FALSE){
  # Set it so numbers up to X are written in full (not scientific notation) - reset at the end
  options(scipen = 1000)
  if(sum(c("Population", "Latitude", "Longitude") %in% colnames(df))<3){
    stop("Input must have at least three columns ; 'Population', 'Latitude' and 'Longitude'")
  }
  if(is.null(size_col)){
    df$size=1
    size='size'
  }else{size=size_col}
  if(is.null(colour)){
    if(is.null(fill_col)){colour='white'}else{colour=plasma(100)}
  }
  if(is.null(fill_col)){
    df$fill=1
    fill='fill'
  }else{fill=fill_col}
  size_range=max(df[,size])-min(df[,size])
  size_breaks=unique(c(min(df[,size]), 
                       min(df[,size])+round(0.125*size_range),
                       min(df[,size])+round(0.25*size_range),
                       min(df[,size])+round(0.5*size_range),
                       #min(df[,size])+round(0.75*size_range),
                       max(df[,size])))
  if(length(size_breaks)>nrow(df)){
    size_breaks=unique(df[,size])[order(unique(df[,size]))]
  }
  fill_range=max(df[,fill])-min(df[,fill])
  fill_breaks=unique(c(min(df[,fill]), 
                       min(df[,fill])+round(0.25*fill_range),
                       min(df[,fill])+round(0.5*fill_range),
                       min(df[,fill])+round(0.75*fill_range),
                       max(df[,fill])))
  if(length(fill_breaks)>nrow(df)){
    fill_breaks=unique(df[,fill])[order(unique(df[,fill]))]
  }
  if(!is.null(fill_col)){
    if(grepl("%", fill_col)){
      #fill_breaks=c(floor(min(df[,fill])), 25, 50, 75, ceiling(max(df[,fill])))
      fill_breaks=c(0, 25, 50, 75, 100)
    }
  }
  # Load maps
  if(is.null(world_map)){
    world_map=map_data("world")
  }
  if(is.null(chimp.ranges)){
    chimp.ranges=readOGR("~/OneDrive - University College London/Projects/Ostridge_PanAf/environmental_data/data/maps/chimp_ranges/chimp_ranges.shp")
  }
  if(is.null(rivers)){
    ## Downloaded version 5.0.0 from https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-rivers-lake-centerlines/
    rivers=readOGR("~/OneDrive - University College London/Projects/Ostridge_PanAf/environmental_data/data/maps/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp")
    rivers=rivers[grepl("Congo|Niger|Ubangi|Ogooué|Tanganyika|Sanaga|Benue|Uelé|Lualaba|Chambeshi|Chinko|Mbomou|Uele|Mbam", rivers$name_en),]
  }
  if(is.null(lakes)){
    ## Download water bodies Zip (last updated Nov 14, 2018) https://datacatalog.worldbank.org/search/dataset/0040797
    lakes=readOGR("~/OneDrive - University College London/Projects/Ostridge_PanAf/environmental_data/data/maps/africawaterbody/Africa_waterbody.shp")
  }
  # Limits
  ## x
  xmin=min(df$Longitude)
  xmax=max(df$Longitude)
  xmin=xmin-(xmax-xmin)*0.05
  xmax=xmax+(xmax-xmin)*0.05
  ## x
  ymin=min(df$Latitude)
  ymax=max(df$Latitude)
  ymin=ymin-(ymax-ymin)*0.15
  ymax=ymax+(ymax-ymin)*0.15
  if(!is.null(ymax_over_ride)){ymax=ymax_over_ride}
  if(!is.null(ymin_over_ride)){ymin=ymin_over_ride}
  # name if plotting forest tree percentage 
  #if(fill_var=="forest_tree_pct"){legend_name="Forest\ntrees (%)"}else{legend_name=env_var}
  # Plot
  r_alpha=1
  ## Scale
  scale_size=1000
  l_r=0.25
  scale_xmin=xmin+l_r*(xmax-xmin)
  # This function is also defined in this script and uses the haversine formula
  scale_xmax=calculate_new_longitude(lon_old=scale_xmin, lat_old=ymin, distance_km=1000, bearing=90)
  scale_xrange=scale_xmax-scale_xmin
  ## Plot
  df=df[order(-df[[size]], df[["Population"]]),]
  map=ggplot(NULL) +
    ## Map
    theme_void() +
    geom_polygon(data=world_map, aes(x = long, y = lat, group = group), fill="navajowhite3", colour = "navajowhite4", alpha=0.25) +
    coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) +
    labs(title=title) +
    theme(plot.title = element_text(hjust = 0.5, size=15),
          legend.key.height= unit(5, 'mm'), legend.key.width= unit(6, 'mm'),
          legend.position = legend.position,
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          legend.direction = "vertical", legend.box = "horizontal", 
          legend.key = element_rect(fill = "transparent", colour=NA),
          legend.box.background = element_rect(fill="transparent",size=0)) +
    ## Subspecies ranges
    ## White underlayer (ensures colours stay the same even when when slightly transparent)
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies %in% c("troglodytes","schweinfurthii", "ellioti", "verus"),], 
                 aes(x = long, y = lat, group = group), fill="white", alpha=1) +
    ## Colours
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="troglodytes",], aes(x = long, y = lat, group = group), fill="green2", alpha=r_alpha) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="schweinfurthii",], aes(x = long, y = lat, group = group), fill="darkorange", alpha=r_alpha) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="ellioti",], aes(x = long, y = lat, group = group), fill="red3", alpha=r_alpha) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="verus",], aes(x = long, y = lat, group = group), fill="blue", alpha=r_alpha) +
    ## Borders
    geom_polygon(data=world_map, aes(x = long, y = lat, group = group), fill="navajowhite3", colour = "navajowhite4", alpha=0) +
    ## Water
    ### White outline
    #geom_polygon(data=lakes, aes(x = long, y = lat, group=group), col="white", alpha=1, size=1)+
    geom_path(data=rivers, aes(x = long, y = lat, group = group), col="white", alpha=1, size=1) +
    ### navajowhite3 outline (on top of white)
    #geom_polygon(data=lakes, aes(x = long, y = lat, group=group), col=alpha("navajowhite3", 0.25), alpha=0.25, size=1)+
    #geom_path(data=rivers, aes(x = long, y = lat, group = group), col="navajowhite3", alpha=0.25, size=1) +
    ### Blue
    geom_polygon(data=lakes, aes(x = long, y = lat, group=group), fill="deepskyblue", alpha=1)+
    geom_path(data=rivers, aes(x = long, y = lat, group = group), col="deepskyblue", alpha=1) +
    ## Population labels
    geom_text_repel(data=df, aes(x=Longitude, y=Latitude, label=Population), 
                    #fontface= "bold",
                    min.segment.length = unit(0, 'lines'),
                    size=text_size,
                    seed=seed,
                    max.overlaps = Inf,
                    bg.color = "white",
                    bg.r = 0.1,
                    fontface = 'bold',
                    force=force,
                    box.padding=unit(box.padding, 'lines')
    ) +
    ### Cover annoying little island (Ascension Island I think) which gets in the way of the legends
    geom_rect(aes(xmin = -15, xmax = -13, ymin = -9, ymax =-6), fill = "white", size = 0)+
    ## Plot points
    geom_point(data=df, aes_string(x='Longitude', y='Latitude', size=paste0('`', size, '`'), fill=paste0('`', fill, '`')), colour='black', pch=21) +
    #scale_colour_gradientn(name=fill, colours = colour, aesthetics = c("fill"), breaks = fill_breaks) +
    scale_size_continuous(range  = c(min_size, max_size), 
                          limits = c(min(df[,size]), max(df[,size])),
                          breaks = size_breaks) +
    guides(colour = guide_legend(order=2),
           size = guide_legend(order=1)) +
    ## Population labels (once more but with no lines so the text overlays the points but lines do not)
    geom_text_repel(data=df, aes(x=Longitude, y=Latitude, label=Population), 
                    #fontface= "bold",
                    min.segment.length = unit(999, 'lines'),
                    size=text_size,
                    seed=seed,
                    max.overlaps = Inf,
                    bg.color = "white",
                    bg.r = 0.1,
                    fontface = 'bold',
                    force=force,
                    box.padding=unit(box.padding, 'lines')
    )
  if(scale){
    # Scales
    map=map+geom_rect(aes(xmin = scale_xmin, xmax = scale_xmax, ymin = ymin, ymax = ymin + 0.03*(ymax-ymin)), fill = "black", colour="black", size = 1)+
      geom_rect(aes(xmin = scale_xmin, xmax = scale_xmin+scale_xrange*0.25, ymin = ymin, ymax = ymin + 0.03*(ymax-ymin)), fill = "white", size = 0)+
      geom_rect(aes(xmin = scale_xmin+scale_xrange*0.5, xmax = scale_xmin+scale_xrange*0.75, ymin = ymin, ymax = ymin + 0.03*(ymax-ymin)), fill = "white", size = 0)+
      
      geom_rect(aes(xmin = scale_xmin+scale_xrange*0.5, xmax = scale_xmin+scale_xrange*0.75, ymin = ymin, ymax = ymin + 0.015*(ymax-ymin)), fill = "black", size = 0)+
      geom_rect(aes(xmin = scale_xmin+scale_xrange*0.75, xmax = scale_xmax, ymin = ymin, ymax = ymin + 0.015*(ymax-ymin)), fill = "white", size = 0)+
      geom_rect(aes(xmin = scale_xmin, xmax = scale_xmin+scale_xrange*0.25, ymin = ymin, ymax = ymin + 0.015*(ymax-ymin)), fill = "black", size = 0)+
      geom_rect(aes(xmin = scale_xmin+scale_xrange*0.25, xmax = scale_xmin+scale_xrange*0.5, ymin = ymin, ymax = ymin + 0.015*(ymax-ymin)), fill = "white", size = 0)+
      
      geom_text(aes(x =  mean(scale_xmax, scale_xmin), y = ymin + 0.06*(ymax-ymin), label = paste0("     ", scale_size, " km")), color = "black", size = 3, fontface = "bold")+
      geom_text(aes(x = scale_xmin, y = ymin + 0.06*(ymax-ymin), label = paste0("0")), color = "black", size = 3, fontface = "bold") +
      # Compass
      #annotation_custom(rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")),
      #                  xmin = scale_xmin, xmax = scale_xmax,
      #                  ymin = ymin, ymax = ymin + 0.4*(ymax-ymin)) +
      geom_polygon(aes(x = c(scale_xmin+scale_xrange*0.5, scale_xmin+scale_xrange*0.5, scale_xmin+scale_xrange*0.65), 
                       y = c(ymin + 0.25*(ymax-ymin), ymin+0.175*(ymax-ymin), ymin + 0.15*(ymax-ymin))), 
                   fill = "black", color = "black")+
      geom_polygon(aes(x = c(scale_xmin+scale_xrange*0.5, scale_xmin+scale_xrange*0.5, scale_xmin+scale_xrange*0.35), 
                       y = c(ymin + 0.25*(ymax-ymin), ymin+0.175*(ymax-ymin), ymin + 0.15*(ymax-ymin))), 
                   fill = "white", color = "black")+
      #geom_segment(aes(x=mean(c(scale_xmax, scale_xmin)), y=ymin+0.15*(ymax-ymin), 
      #                 xend=mean(c(scale_xmax, scale_xmin)), yend=ymin + 0.25*(ymax-ymin)),
      #             arrow = arrow(length = unit(0.4, "cm")), size=1) +
      geom_text(aes(x =mean(c(scale_xmax, scale_xmin)), y = ymin+0.125*(ymax-ymin), label = "N"), color = "black", size = 5, fontface = "bold")
  }
  #if(!is.null(fill_col)){
    if(log_fill){
      map=map+scale_colour_gradientn(name=fill, colours = colour, aesthetics = c("fill"), trans = "log10", limits = c(min(fill_breaks), max(fill_breaks)))
    }else{
      map=map+scale_colour_gradientn(name=fill, colours = colour, aesthetics = c("fill"), breaks = fill_breaks, limits = c(min(fill_breaks), max(fill_breaks)))
    }
  #}else{
  if(is.null(fill_col)){
    map=map+guides(fill = F)
  }
  if(is.null(size_col)){
    map=map+guides(size = F)
  }
  print(map)
  # Reset
  options("scipen") 
}

plot_fpr_stats_main=function(xtx_fpr, null, fprs=c(0.005, 0.001, 0.0005), thresh_stat='`XtXst.med`'){
  subsps=names(xtx_fpr)
  fdr.df=NULL
  cat("Coverage Bin Stats\n")
  for(subsp in subsps){
    cat(subsp, "\n")
    for(fpr in fprs){
      bin_df=data.frame(fpr=numeric(), bin=numeric(), COVARIABLE_name=character(), n_SNPs=numeric(), dataset=character(), thresh=numeric())
      # Get total number of SNPs
      tot_n_snps=nrow(unique(xtx_fpr[[subsp]][,c('chr', 'pos')]))
      # Select candidates
      cands=unique(xtx_fpr[[subsp]][xtx_fpr[[subsp]]$fpr<=fpr, ])
      # Get number of candidate SNPs
      cand_n_snps=nrow(cands)
      # FDR = number of false positives / total number of positives.
      ## number of false positives = FPR x total number SNPs.
      ## total number of positives = number of candidates.
      fdr=(fpr*tot_n_snps)/cand_n_snps
      ## FDR > 1 doesn't really make sense so I round it down to 1
      fdr=min(1, fdr)
      # Make data frame for plotting number of candidates 
      row=data.frame('Subspecies'=subsp, 
                     'FPR'=fpr,
                     'Total_SNPs'= tot_n_snps,
                     'Candidate_SNPs'=cand_n_snps,
                     'Expected_SNPs'=fpr*tot_n_snps,
                     'Threshold'= min(cands[[thresh_stat]]),
                     'FDR'=fdr,
                     'True_Positives'=cand_n_snps*(1-fdr),
                     'False_Positives'=cand_n_snps*fdr)
      if(!is.null(fdr.df)){
        fdr.df=rbind(fdr.df, row)
      }else{
        fdr.df=row
      }
    }
  }
  # Plot number of candidates and expected number
  cat("Number of Candidates\n")
  fdr.df=fdr.df[, c('Subspecies', 'FPR', 'Candidate_SNPs', 'Expected_SNPs')]
  # For each fpr value
  # Plot number positives and null expectation for a specific FPR
  ## Rename so plot looks nicer
  fdr.df[fdr.df$Subspecies=='all', 'Subspecies']="All Subspecies"
  fdr.df[fdr.df$Subspecies=='c', 'Subspecies']="Central"
  fdr.df[fdr.df$Subspecies=='e', 'Subspecies']="Eastern"
  fdr.df[fdr.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
  fdr.df[fdr.df$Subspecies=='n', 'Subspecies']="Nigeria-Cameroon"
  fdr.df[fdr.df$Subspecies=='w', 'Subspecies']="Western"
  #fdr.df$Subspecies=factor(fdr.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
  ## Plot
  ### Format file for input
  fdr.df$FPR=factor(paste0(100*fdr.df$FPR, "%"), levels=c('0.5%', '0.1%', '0.05%'))
  bar_plot=ggplot(fdr.df, aes(x=FPR, y=Candidate_SNPs, fill=Subspecies)) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.text.x = element_text(face="bold"),
          strip.background = element_rect(fill="white"),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    geom_bar(position="dodge", stat="identity", col='black', alpha=1) +
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs), size = 1.4, col="white", width=0.85) + 
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs, col="Null\nexpectation"), size = 1, width=0.8) +
    facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies)), scales = "free") +
    #facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies))) +
    labs(x="Estimated FPR threshold", y="Number of SNPs") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_discrete_manual(name='', values = c("Null\nexpectation"='violetred'), 
                          aesthetics = c("col")) +
    scale_discrete_manual(values = c("Central"='green3', "Eastern"='darkorange', "Central-Eastern"='brown', "Nigeria-Cameroon"='red', "Western"='blue', "All Subspecies"="lightblue4"), 
                          aesthetics = c("fill"), guide="none") 
  print(bar_plot)
  return(fdr.df)
}

plot_fpr_stats_main.v2=function(xtx_fpr, null, fprs=c(0.005, 0.001, 0.0005), thresh_stat='`XtXst.med`'){
  subsps=names(xtx_fpr)
  fdr.df=NULL
  cat("Coverage Bin Stats\n")
  for(subsp in subsps){
    cat(subsp, "\n")
    for(fpr in fprs){
      bin_df=data.frame(fpr=numeric(), bin=numeric(), COVARIABLE_name=character(), n_SNPs=numeric(), dataset=character(), thresh=numeric())
      # Get total number of SNPs
      tot_n_snps=nrow(unique(xtx_fpr[[subsp]][,c('chr', 'pos')]))
      # Select candidates
      cands=unique(xtx_fpr[[subsp]][xtx_fpr[[subsp]]$fpr<=fpr, ])
      # Get number of candidate SNPs
      cand_n_snps=nrow(cands)
      # FDR = number of false positives / total number of positives.
      ## number of false positives = FPR x total number SNPs.
      ## total number of positives = number of candidates.
      fdr=(fpr*tot_n_snps)/cand_n_snps
      ## FDR > 1 doesn't really make sense so I round it down to 1
      fdr=min(1, fdr)
      # Make data frame for plotting number of candidates 
      row=data.frame('Subspecies'=subsp, 
                     'FPR'=fpr,
                     'Total_SNPs'= tot_n_snps,
                     'Candidate_SNPs'=cand_n_snps,
                     'Expected_SNPs'=fpr*tot_n_snps,
                     'Threshold'= min(cands[[thresh_stat]]),
                     'FDR'=fdr,
                     'True_Positives'=cand_n_snps*(1-fdr),
                     'False_Positives'=cand_n_snps*fdr)
      if(!is.null(fdr.df)){
        fdr.df=rbind(fdr.df, row)
      }else{
        fdr.df=row
      }
    }
  }
  # Plot number of candidates and expected number
  cat("Number of Candidates\n")
  fdr.df=fdr.df[, c('Subspecies', 'FPR', 'Candidate_SNPs', 'Expected_SNPs')]
  # For each fpr value
  # Plot number positives and null expectation for a specific FPR
  ## Rename so plot looks nicer
  fdr.df[fdr.df$Subspecies=='all', 'Subspecies']="All Subspecies"
  fdr.df[fdr.df$Subspecies=='c', 'Subspecies']="Central"
  fdr.df[fdr.df$Subspecies=='e', 'Subspecies']="Eastern"
  fdr.df[fdr.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
  fdr.df[fdr.df$Subspecies=='n', 'Subspecies']="Nigeria-Cameroon"
  fdr.df[fdr.df$Subspecies=='w', 'Subspecies']="Western"
  #fdr.df$Subspecies=factor(fdr.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
  ## Plot
  ### Format file for input
  fdr.df$FPR=factor(paste0(100*fdr.df$FPR, "%"), levels=c('0.5%', '0.1%', '0.05%'))
  bar_plot=ggplot(fdr.df, aes(x=FPR, y=Candidate_SNPs, fill=Subspecies)) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.text.x = element_text(face="bold"),
          strip.background = element_rect(fill="white"),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    geom_bar(position="dodge", stat="identity", col='black', alpha=1) +
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs), size = 2, col="white", width=0.85) + 
    #geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs, col="Null\nexpectation"), size = 1, width=0.8) +
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs, col=Subspecies), size = 1, width=0.8, linetype="solid") +
    facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies)), scales = "free") +
    #facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies))) +
    labs(x="Estimated FPR threshold", y="Number of SNPs") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    #cale_discrete_manual(name='', values = c("Null\nexpectation"='violetred'), 
    #                      aesthetics = c("col")) +
    scale_discrete_manual(values = c("Central"='green3', "Eastern"='darkorange', "Central-Eastern"='maroon', "Nigeria-Cameroon"='red', "Western"='blue', "All Subspecies"="lightblue4"), 
                          aesthetics = c("col", "fill"), guide="none") 
  print(bar_plot)
  return(fdr.df)
}

plot_fpr_stats_main.v3=function(xtx_fpr, null, fprs=c(0.005, 0.001, 0.0005), thresh_stat='XtXst.med', title=NULL){
  if(is.null(title)){
    title=paste0("")
  }
  #xtx_fpr=unique(xtx_fpr[colnames(xtx_fpr)!='gene',])
  subsps=names(xtx_fpr)
  fdr.df=NULL
  cat("Coverage Bin Stats\n")
  for(subsp in subsps){
    cat(subsp, "\n")
    for(fpr in fprs){
      bin_df=data.frame(fpr=numeric(), bin=numeric(), COVARIABLE_name=character(), n_SNPs=numeric(), dataset=character(), thresh=numeric())
      # Get total number of SNPs
      tot_n_snps=nrow(unique(xtx_fpr[[subsp]][,c('chr', 'pos')]))
      # Select candidates
      cands=unique(xtx_fpr[[subsp]][xtx_fpr[[subsp]]$fpr<=fpr, c('chr', 'pos', thresh_stat)])
      # Get number of candidate SNPs
      cand_n_snps=nrow(cands)
      # FDR = number of false positives / total number of positives.
      ## number of false positives = FPR x total number SNPs.
      ## total number of positives = number of candidates.
      fdr=(fpr*tot_n_snps)/cand_n_snps
      ## FDR > 1 doesn't really make sense so I round it down to 1
      fdr=min(1, fdr)
      # Make data frame for plotting number of candidates 
      row=data.frame('Subspecies'=subsp, 
                     'FPR'=fpr,
                     'Total_SNPs'= tot_n_snps,
                     'Candidate_SNPs'=cand_n_snps,
                     'Expected_SNPs'=fpr*tot_n_snps,
                     'Threshold'= min(cands[[thresh_stat]]),
                     'FDR'=fdr,
                     'True_Positives'=cand_n_snps*(1-fdr),
                     'False_Positives'=cand_n_snps*fdr)
      if(!is.null(fdr.df)){
        fdr.df=rbind(fdr.df, row)
      }else{
        fdr.df=row
      }
    }
  }
  # Plot number of candidates and expected number
  cat("Number of Candidates\n")
  fdr.df=fdr.df[, c('Subspecies', 'FPR', 'Candidate_SNPs', 'Expected_SNPs')]
  # For each fpr value
  # Plot number positives and null expectation for a specific FPR
  ## Rename so plot looks nicer
  fdr.df[fdr.df$Subspecies=='all', 'Subspecies']="All"
  fdr.df[fdr.df$Subspecies=='c', 'Subspecies']="Central"
  fdr.df[fdr.df$Subspecies=='e', 'Subspecies']="Eastern"
  fdr.df[fdr.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
  fdr.df[fdr.df$Subspecies=='n', 'Subspecies']="Nigeria-Cameroon"
  fdr.df[fdr.df$Subspecies=='w', 'Subspecies']="Western"
  #fdr.df$Subspecies=factor(fdr.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
  ## Plot
  ### Format file for input
  #fdr.df$FPR=factor(paste0(100*fdr.df$FPR, "%"), levels=c('0.5%', '0.1%', '0.05%'))
  fdr.df$FPR=factor(paste0(100*fdr.df$FPR, "%"), levels=paste0(fprs*100,'%'))
  bar_plot=ggplot(fdr.df, aes(x=FPR, y=Candidate_SNPs, fill=Subspecies)) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          #strip.text.x = element_text(face="bold"),
          strip.text = element_text(face = "italic"),
          strip.background = element_rect(fill="white"),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    geom_bar(position="dodge", stat="identity", col='black', alpha=1) +
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs), size = 2, col="black", width=0.85) + 
    #geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs, col="Null\nexpectation"), size = 1, width=0.8) +
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs, col=Subspecies), size = 0.5, width=0.75, linetype="solid", col="white") +
    facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies)), scales = "free") +
    #facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies))) +
    labs(title=title,
         x="Estimated FPR threshold", y="Number of SNPs") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    #cale_discrete_manual(name='', values = c("Null\nexpectation"='violetred'), 
    #                      aesthetics = c("col")) +
    scale_discrete_manual(values = c("Central"='green3', "Eastern"='darkorange', 
                                     #"Central-Eastern"="khaki4", 
                                     "Central-Eastern"="gold", 
                                     "Nigeria-Cameroon"='red', "Western"='blue', "All"="lightblue4"), 
                          aesthetics = c("col", "fill"), guide="none") 
  print(bar_plot)
  return(fdr.df)
}

plot_fpr_stats_main.v4=function(xtx_fpr, null, fprs=c(0.005, 0.001, 0.0005), thresh_stat='XtXst.med', title=NULL, subtitle=NULL, verbose=TRUE, plot_ratio=TRUE, tick_int=NULL){
  if(is.null(title)){
    title=paste0("")
  }
  #xtx_fpr=unique(xtx_fpr[colnames(xtx_fpr)!='gene',])
  subsps=names(xtx_fpr)
  fdr.df=NULL
  if(verbose){cat("Coverage Bin Stats\n")}
  for(subsp in subsps){
    if(verbose){cat(subsp, "\n")}
    for(fpr in fprs){
      bin_df=data.frame(fpr=numeric(), bin=numeric(), COVARIABLE_name=character(), n_SNPs=numeric(), dataset=character(), thresh=numeric())
      # Get total number of SNPs
      tot_n_snps=nrow(unique(xtx_fpr[[subsp]][,c('chr', 'pos')]))
      # Select candidates
      cands=unique(xtx_fpr[[subsp]][xtx_fpr[[subsp]]$fpr<=fpr, c('chr', 'pos', thresh_stat)])
      # Get number of candidate SNPs
      cand_n_snps=nrow(cands)
      # FDR = number of false positives / total number of positives.
      ## number of false positives = FPR x total number SNPs.
      ## total number of positives = number of candidates.
      fdr=(fpr*tot_n_snps)/cand_n_snps
      ## FDR > 1 doesn't really make sense so I round it down to 1
      fdr=min(1, fdr)
      # Make data frame for plotting number of candidates 
      row=data.frame('Subspecies'=subsp, 
                     'FPR'=fpr,
                     'Total_SNPs'= tot_n_snps,
                     'Candidate_SNPs'=cand_n_snps,
                     'Expected_SNPs'=fpr*tot_n_snps,
                     'Ratio'=round(cand_n_snps/(fpr*tot_n_snps),2),
                     'Threshold'= min(cands[[thresh_stat]]),
                     'FDR'=fdr,
                     'True_Positives'=cand_n_snps*(1-fdr),
                     'False_Positives'=cand_n_snps*fdr)
      if(!is.null(fdr.df)){
        fdr.df=rbind(fdr.df, row)
      }else{
        fdr.df=row
      }
    }
  }
  # Plot number of candidates and expected number
  if(verbose){cat("Number of Candidates\n")}
  fdr.df=fdr.df[, c('Subspecies', 'FPR', 'Candidate_SNPs', 'Expected_SNPs', 'Ratio')]
  # For each fpr value
  # Plot number positives and null expectation for a specific FPR
  ## Rename so plot looks nicer
  fdr.df[fdr.df$Subspecies=='all', 'Subspecies']="All"
  fdr.df[fdr.df$Subspecies=='c', 'Subspecies']="Central"
  fdr.df[fdr.df$Subspecies=='e', 'Subspecies']="Eastern"
  fdr.df[fdr.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
  fdr.df[fdr.df$Subspecies=='n', 'Subspecies']="Nigeria-Cameroon"
  fdr.df[fdr.df$Subspecies=='w', 'Subspecies']="Western"
  #fdr.df$Subspecies=factor(fdr.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
  ## Plot
  ### Format file for input
  #fdr.df$FPR=factor(paste0(100*fdr.df$FPR, "%"), levels=c('0.5%', '0.1%', '0.05%'))
  fdr.df$FPR=factor(paste0(100*fdr.df$FPR, "%"), levels=paste0(fprs*100,'%'))
  tmp=fdr.df %>%
    group_by(Subspecies) %>%
    mutate(y_max = pmax(Candidate_SNPs, Expected_SNPs)) %>%
    slice_max(y_max) %>%
    dplyr::select(Subspecies, y_max)
  tmp$y_max=1.075*tmp$y_max
  fdr.df=merge(fdr.df, tmp, by='Subspecies')
  
  bar_plot=ggplot(fdr.df, aes(x=FPR, y=Candidate_SNPs, fill=Subspecies)) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          #strip.text.x = element_text(face="bold"),
          strip.text = element_text(face = "italic"),
          strip.background = element_rect(fill="white"),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    geom_bar(position="dodge", stat="identity", col='black', alpha=1) +
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs), size = 2, col="black", width=0.85) + 
    #geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs, col="Null\nexpectation"), size = 1, width=0.8) +
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs, col=Subspecies), size = 0.5, width=0.75, linetype="solid", col="white") +
    facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies)), scales = "free") +
    #facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies))) +
    labs(title=title, subtitle=subtitle,
         x="Estimated FPR threshold", y="Number of SNPs") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    #scale_discrete_manual(name='', values = c("Null\nexpectation"='violetred'), 
    #                      aesthetics = c("col")) +
    scale_discrete_manual(values = c("Central"='green3', "Eastern"='darkorange', 
                                     #"Central-Eastern"="khaki4", 
                                     "Central-Eastern"="gold", 
                                     "Nigeria-Cameroon"='red', "Western"='blue', "All"="lightblue4"), 
                          aesthetics = c("col", "fill"), guide="none") 
  if(!is.null(tick_int)){
    bar_plot=bar_plot+
      # Set consistent y axis ticks
      scale_y_continuous(breaks = seq(0, 100000, by = tick_int)) +
      theme(panel.grid.minor = element_blank())
  }
  if(plot_ratio){
    bar_plot=bar_plot+
      geom_text(aes(y=y_max, label=Ratio), col="black", fontface="bold", size=3)
  }
  print(bar_plot)
  return(fdr.df)
}
