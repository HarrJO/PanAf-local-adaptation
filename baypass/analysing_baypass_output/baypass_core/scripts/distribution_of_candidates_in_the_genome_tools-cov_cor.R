library(ggplot2)
library(wesanderson)

# Plot
manhattan_plot_v4=function(xtx, title="Manhattan Plot", stat="XtXst.med", chr.lengths=NULL, gtf=NULL, focal.snps=NULL, p.thresh=c(0.005, 0.001, 0.0005), lines=FALSE, ymin=NULL){
  xtx=xtx[order(xtx$chr),]
  xtx$chr=as.factor(xtx$chr)
  chrs=unique(xtx$chr)
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
  ## Add cumulative position to xtx data
  xtx=merge(xtx, cum.chr.start, by='chr')
  xtx$cumpos=as.numeric(xtx$pos + xtx$cum.chr.start)
  # Prepare table for plotting the positions of genes if a gtf file is provided
  if(!is.null(gtf)){
    ## Keep only selected chrs
    gtf=gtf[gtf$V1 %in% chrs,]
    ## Get the maximum and minimum positions of the gene (as the gtf may contain CDS and exon information instead of whole gene coordinates)
    ### V10 is the gene name column, V1 is chr, V4 is start and V5 is end
    gene.bed=gtf %>% group_by(V10) %>% summarize(chr=unique(V1), start=min(V4), end=max(V5))
    colnames(gene.bed)[1]='gene'
    ## Add cumulative position of chr
    gene.bed=merge(gene.bed, cum.chr.start, by='chr')
    ## Make cumulative start and ends for gene
    gene.bed$cumstart=gene.bed$start + gene.bed$cum.chr.start
    gene.bed$cumend=gene.bed$end + gene.bed$cum.chr.start
  }
  # Manhatten Plot
  ## Get thresholds
  stat.thresh=c()
  for(p in p.thresh){
    stat.thresh=c(stat.thresh, min(xtx[xtx$empirical_p_med_cov<=p, stat]))
  }
  ## Prepare x axis
  ### Midpoints
  chr.mid=data.frame(chrs=chrs, center=(cumsum+cum.chr.start$cum.chr.start)/2)
  ### Min and Max
  chr.max=cumsum
  chr.min=cum.chr.start$cum.chr.start
  ## Plotting
  ### Set plotting colours
  thresh.col=wes_palette("Zissou1", n = length(p.thresh), type = "continuous")
  chr.panel.col=rep(c("grey", "white"), 20)
  chr.point.col=rep(c("cyan4", "sienna1"), 20)
  #### Crop to the number of chrs
  chr.panel.col=chr.panel.col[1:length(chrs)]
  chr.point.col=chr.point.col[1:length(chrs)]
  ### Plot
  man.plot=ggplot(NULL) +
    #### General
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour = chr.panel.col, face="bold", size=7)) +
    labs(title=paste0(title), x="Chromosome", y=paste0(stat)) +
    #### Axes
    scale_x_continuous(label = chr.mid$chr, breaks= chr.mid$center ) +
    scale_y_continuous(expand = c(0, 0) ) +
    geom_rect(aes(xmin=chr.min, xmax=chr.max, ymin=-Inf, max=Inf), fill=chr.panel.col, alpha=0.5) +
    #### Add thresholds
    geom_hline(yintercept=stat.thresh, linetype="dashed", colour=thresh.col) +
    geom_label(data=NULL,aes(max(xtx$cumpos)*1.07, stat.thresh, label = paste0(100*p.thresh,"%")), colour=thresh.col, fontface = "bold")
  #### Add gene positions?
  if(!is.null(gtf)){
    man.plot=man.plot+geom_rect(aes(xmin=gene.bed$cumstart, xmax=gene.bed$cumend, ymin=-Inf, max=Inf), fill='red', alpha=max(0.25,1/nrow(gene.bed)))
  }
  #### Add SNP stats
  cumpos="cumpos"
  chr="chr"
  ##### Add dots or lines connecting points?
  if(lines){
    man.plot=man.plot+geom_line(data=xtx, aes_string(x=cumpos, y=stat, col=chr), size=0.1)
  }else{
    man.plot=man.plot+
      geom_point(data=xtx, aes_string(x=cumpos, y=stat, col=chr), alpha=0.5, size=0.5)
  }
  man.plot=man.plot + scale_color_manual(values = chr.point.col)
  ### Highlight SNPs?
  if(!is.null(focal.snps)){
    ### Get the results for focal SNPs (doing a left merge (i.e. all.x=TRUE) we only keep results for our focal SNPs)
    focal.snps.xtx=merge(focal.snps[,c('chr','pos')], xtx, by=c('chr', 'pos'), all.x=TRUE)
    ### Add to plot
    man.plot=man.plot+geom_point(data=focal.snps.xtx, aes_string(x=cumpos, y=stat), colour="black", alpha=1, size=0.5)
  }
  #### Add gene names if not many
  if(!is.null(gtf)){
    if(nrow(gene.bed)<20){
      ##### For some reason you need to assign this to an object first then put it into geom_label()
      max.y=layer_scales(man.plot)$y$range$range[2]
      min.y=layer_scales(man.plot)$y$range$range[1]
      range.y=max.y-min.y
      ##### runif randomly selects a y value in the range  to minimise overlap
      man.plot=man.plot+geom_text(data=gene.bed, aes(x=cumstart, y =runif(nrow(gene.bed),min=min.y+range.y*0.25, max.y-range.y*0.1), label=gene), colour="red", aplha=0.5, size=3)
    }
  }
  #### ylim?
  if(!is.null(ymin)){
    man.plot=man.plot+scale_y_continuous(limits = c(ymin, NA), expand = c(0, 0) )
  }
  print(man.plot)
}


# PLOT SNP DENISTY
SNP_distribution=function(xtx.ann, tails=c(1, 0.005, 0.001, 0.0005), chrs=1:22, chr.lengths=NULL, gtf=NULL){
  xtx.ann=as.data.frame(xtx.ann)
  xtx.ann=xtx.ann[!is.na(xtx.ann$gene),]
  xtx=unique(xtx.ann[colnames(xtx.ann)!='gene'])
  # Make sure chrs are in ascending order (this ensures all other data frames are in the same order)
  chrs=chrs[order(chrs)]
  # Only include selected chromosomes
  xtx.ann=xtx.ann[xtx.ann$chr %in% chrs,]
  xtx=xtx[xtx$chr %in% chrs,]
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
  ## Add cumulative position to xtx data
  xtx=merge(xtx, cum.chr.start, by='chr')
  xtx$cumpos=xtx$pos + xtx$cum.chr.start
  # Prepare table for plotting the positions of genes if a gtf file is provided
  if(!is.null(gtf)){
    ## Keep only selected chrs
    gtf=gtf[gtf$V1 %in% chrs,]
    ## Get the maximum and minimum positions of the gene (as the gtf may contain CDS and exon information instead of whole gene coordinates)
    ### V10 is the gene name column, V1 is chr, V4 is start and V5 is end
    gene.bed=gtf %>% group_by(V10) %>% summarize(chr=unique(V1), start=min(V4), end=max(V5))
    colnames(gene.bed)[1]='gene'
    ## Add cumulative position of chr
    gene.bed=merge(gene.bed, cum.chr.start, by='chr')
    ## Make cumulative start and ends for gene
    gene.bed$cumstart=gene.bed$start + gene.bed$cum.chr.start
    gene.bed$cumend=gene.bed$end + gene.bed$cum.chr.start
  }
  # SNPs and genes
  for(tail in tails){
    ## Select tail SNPs
    #xtx.ann.tail=xtx.ann[xtx.ann$empirical_p_med_cov<tail,]
    xtx.ann.tail=xtx.ann[xtx.ann$empirical_p_med_cov<=tail,]
    #xtx.tail=xtx[xtx$empirical_p_med_cov<tail,]
    xtx.tail=xtx[xtx$empirical_p_med_cov<=tail,]
    ## Number of SNPs
    cat(tail*100, "% tail: \n", length(unique(xtx.tail$MRK)), " 'genic' SNPs\n ")
    ## Number of genes
    cat(length(unique(xtx.ann.tail$gene)), " genes\n ")
    ## SNPs per gene
    cat("Mean SNPs per gene: ", length(unique(xtx.ann.tail$MRK))/length(unique(xtx.ann.tail$gene)), " genes\n")
    snp.per.gene_hist.data=table(xtx.ann.tail$gene)
    snp.per.gene_hist=ggplot(NULL) +
      theme_bw() +
      geom_density(aes(x=snp.per.gene_hist.data), fill='black', alpha=0.5) +
      labs(title=paste0(tail*100, "% tail: number of SNPs per gene"), x="Number of SNPs", y="Frequency of Genes") +
      theme(plot.title = element_text(hjust = 0.5))
    print(snp.per.gene_hist)
    ## Genes per SNP 
    gene.per.snp_hist.data=table(xtx.ann.tail$MRK)
    gene.per.snp_hist=ggplot(NULL) +
      theme_bw() +
      geom_histogram(aes(x=gene.per.snp_hist.data)) +
      labs(title=paste0(tail*100, "% tail: number of genes per SNP"), x="Number of Genes", y="Frequency of SNPs") +
      theme(plot.title = element_text(hjust = 0.5))
    print(gene.per.snp_hist)
    # SNP density plot
    ## Prepare x axis
    ### Midpoints
    chr.mid=data.frame(chrs=chrs, center=(cumsum+cum.chr.start$cum.chr.start)/2)
    ### Min and Max
    chr.max=cumsum
    chr.min=cum.chr.start$cum.chr.start
    ## Plotting
    ### Set plotting colours
    chr.col=rep(c("turquoise3", "white"), 20)
    #### Crop to the number of chrs
    chr.col=chr.col[1:length(chrs)]
    ### Plot
    dens.plot=ggplot(NULL) +
      #### General
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(colour = chr.col, face="bold", size=7)) +
      labs(title=paste0(tail*100, "% tail: density of SNPs over genome"), x="Chromosome", y="SNP Denisty") +
      #### Axes
      scale_x_continuous(label = chr.mid$chr, breaks= chr.mid$center ) +
      scale_y_continuous(expand = c(0, 0) ) +
      geom_rect(aes(xmin=chr.min, xmax=chr.max, ymin=-Inf, max=Inf), fill=chr.col, alpha=0.5)
    #### Add gene positions?
    if(!is.null(gtf)){
      dens.plot=dens.plot+geom_rect(aes(xmin=gene.bed$cumstart, xmax=gene.bed$cumend, ymin=-Inf, max=Inf), fill='red', alpha=max(0.25,1/nrow(gene.bed)))
    }
    #### Add SNP densities
    dens.plot=dens.plot+geom_density(aes(x=xtx.tail$cumpos), col="black", fill="black", alpha=1, bw=1) # bw >10000000 begins to find or general trends
    #### Add gene names if not many
    if(!is.null(gtf)){
      if(nrow(gene.bed)<20){
        ##### For some reason you need to assign this to an object first then put it into geom_label()
        max.y=layer_scales(dens.plot)$y$range$range[2]
        min.y=layer_scales(dens.plot)$y$range$range[1]
        range.y=max.y-min.y
        ##### runif randomly selects a y value in the range  to minimise overlap
        dens.plot=dens.plot+geom_text(data=gene.bed, aes(x=cumstart, y =runif(nrow(gene.bed),min=min.y+range.y*0.25, max.y-range.y*0.1), label=gene), colour="red", aplha=0.5, size=3)
      }
    }
    print(dens.plot)
  }
}

SNP_distribution.v2=function(input, chr.lengths=NULL, gtf=NULL){
  input=input[c('chr', 'pos')]
  chrs=unique(input$chr)
  # Make sure chrs are in ascending order (this ensures all other data frames are in the same order)
  chrs=chrs[order(chrs)]
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
  ## Add cumulative position to data
  input=merge(input, cum.chr.start, by='chr')
  input$cumpos=input$pos + input$cum.chr.start
  # Prepare table for plotting the positions of genes if a gtf file is provided
  if(!is.null(gtf)){
    ## Keep only selected chrs
    gtf=gtf[gtf$V1 %in% chrs,]
    ## Get the maximum and minimum positions of the gene (as the gtf may contain CDS and exon information instead of whole gene coordinates)
    ### V10 is the gene name column, V1 is chr, V4 is start and V5 is end
    gene.bed=gtf %>% group_by(V10) %>% summarize(chr=unique(V1), start=min(V4), end=max(V5))
    colnames(gene.bed)[1]='gene'
    ## Add cumulative position of chr
    gene.bed=merge(gene.bed, cum.chr.start, by='chr')
    ## Make cumulative start and ends for gene
    gene.bed$cumstart=gene.bed$start + gene.bed$cum.chr.start
    gene.bed$cumend=gene.bed$end + gene.bed$cum.chr.start
  }
  # SNP density plot
  ## Prepare x axis
  ### Midpoints
  chr.mid=data.frame(chrs=chrs, center=(cumsum+cum.chr.start$cum.chr.start)/2)
  ### Min and Max
  chr.max=cumsum
  chr.min=cum.chr.start$cum.chr.start
  ## Plotting
  ### Set plotting colours
  chr.col=rep(c("turquoise3", "white"), 20)
  #### Crop to the number of chrs
  chr.col=chr.col[1:length(chrs)]
  ### Plot
  dens.plot=ggplot(NULL) +
    #### General
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(colour = chr.col, face="bold", size=7)) +
    labs(title=paste0("Density of SNPs over genome"), x="Chromosome", y="SNP Denisty") +
    #### Axes
    scale_x_continuous(label = chr.mid$chr, breaks= chr.mid$center ) +
    scale_y_continuous(expand = c(0, 0) ) +
    geom_rect(aes(xmin=chr.min, xmax=chr.max, ymin=-Inf, max=Inf), fill=chr.col, alpha=0.5)
  #### Add gene positions?
  if(!is.null(gtf)){
    dens.plot=dens.plot+geom_rect(aes(xmin=gene.bed$cumstart, xmax=gene.bed$cumend, ymin=-Inf, max=Inf), fill='red', alpha=max(0.25,1/nrow(gene.bed)))
  }
  #### Add SNP densities
  dens.plot=dens.plot+geom_density(aes(x=input$cumpos), col="black", fill="black", alpha=1, bw=1) # bw >10000000 begins to find or general trends
  #### Add gene names if not many
  if(!is.null(gtf)){
    if(nrow(gene.bed)<20){
      ##### For some reason you need to assign this to an object first then put it into geom_label()
      max.y=layer_scales(dens.plot)$y$range$range[2]
      min.y=layer_scales(dens.plot)$y$range$range[1]
      range.y=max.y-min.y
      ##### runif randomly selects a y value in the range  to minimise overlap
      dens.plot=dens.plot+geom_text(data=gene.bed, aes(x=cumstart, y =runif(nrow(gene.bed),min=min.y+range.y*0.25, max.y-range.y*0.1), label=gene), colour="red", aplha=0.5, size=3)
    }
    if(length(chrs)==1){
      dens.plot=dens.plot+xlim(min(input$pos), max(input$pos))
    }
    return(dens.plot)
  }
}



snps.per.gene.across.p=function(xtx.ann, empirical_ps=seq(0, 0.02, 0.0005)){
  # Only include SNPs in genes
  xtx.ann=xtx.ann[!is.na(xtx.ann$gene),]
  # Create empty output dataframe
  out.df=data.frame(tail=numeric(), SNPs.per.genes_data=numeric(), SNPs.per.genes_null=numeric())
  # For each empirical p...
  for(p in empirical_ps){
    ## Data
    ### Select tail SNPs
    #xtx.ann.tail=xtx.ann[xtx.ann$empirical_p_med_cov<p,]
    xtx.ann.tail=xtx.ann[xtx.ann$empirical_p_med_cov<p,]
    ### Get the number of tail SNPs
    n.snps.tail=length(unique(xtx.ann.tail$MRK))
    ### Get the number of genes which overlap tail SNPs
    n.genes.tail=length(unique(xtx.ann.tail$gene))
    ## Null
    n.genes.null=c()
    ### Do the following 10 tmes and use the mean
    for(i in 1:10){
      #### Randomly sample the same number of SNPs in the tail from all sites
      rand.MRK=unique(xtx.ann$MRK)[sample(unique(xtx.ann$MRK), n.snps.tail)]
      xtx.ann.rand=xtx.ann[xtx.ann$MRK %in% rand.MRK,]
      #### Get number of genes which overlap with randomly sampled SNPs
      n.genes.null[i]=length(unique(xtx.ann.rand$gene))
    }
    #### Get the mean
    n.genes.null=mean(n.genes.null)
    ## Add to output
    out.df=rbind(out.df, list(tail=p, SNPs.per.genes_data=n.snps.tail/n.genes.tail, SNPs.per.genes_null=n.snps.tail/n.genes.null))
  }
  # Plot
  plot=ggplot(out.df) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_smooth(aes(x=tail, y=SNPs.per.genes_data, col="Data"), se=F) +
    geom_point(aes(x=tail, y=SNPs.per.genes_data, col="Data")) +
    geom_smooth(aes(x=tail, y=SNPs.per.genes_null, col="Null"), se=F) +
    geom_point(aes(x=tail, y=SNPs.per.genes_null, col="Null")) +
    scale_discrete_manual(name='', values = c('Data' = 'red', 'Null' = 'blue'), aesthetics='colour') +
    labs(title="Genes per SNP across emirical p-value thresholds\ncompared to random null",x="Empirical p-value threshold", y = "# SNPs/# Genes")
  print(plot)
}

distance_between_cand_snps_1=function(xtx, p, bw=1){
  # Make sure there is a row per SNP
  xtx=unique(xtx[colnames(xtx)!='gene'])
  ## Data
  ### Select tail SNPs
  xtx.tail=xtx[xtx$empirical_p_med_cov<=p,]
  ## For each chr
  chrs=unique(xtx.tail$chr)
  tail_dists=c()
  tail_dists_min=c()
  for(chr in chrs[order(chrs)]){
    # Distance of tails SNPs
    tail_dists_chr=dist(xtx.tail[xtx.tail$chr==chr,'pos'])
    tail_dists=c(tail_dists, tail_dists_chr)
    ## Get minimum distances
    tail_dists_chr_m=as.matrix(tail_dists_chr)
    tail_dists_chr_m[tail_dists_chr_m==0]=NA
    tail_dists_chr_min=apply(tail_dists_chr_m, 1, FUN = min, na.rm = TRUE)
    tail_dists_min=c(tail_dists_min, tail_dists_chr_min)
  }
  # Null SNPs
  null_dists=list(1,2,3,4,5)
  null_dists_min=list(1,2,3,4,5)
  for(i in 1:length(null_dists)){
    for(chr in chrs[order(chrs)]){
      ### Randomly sample the same number of SNPs in the tail from all sites
      xtx.null=xtx[sample(nrow(xtx), nrow(xtx.tail)),]
      # Distance of null SNPs
      null_dists_chr=dist(xtx.null[xtx.null$chr==chr,'pos'])
      null_dists[[i]]=c(null_dists[[i]], null_dists_chr)
      ## Get minimum distances
      null_dists_chr_m=as.matrix(null_dists_chr)
      null_dists_chr_m[null_dists_chr_m==0]=NA
      null_dists_chr_min=apply(null_dists_chr_m, 1, FUN = min, na.rm = TRUE)
      null_dists_min[[i]]=c(null_dists_min[[i]], null_dists_chr_min)
    }
  }
  # Plot
  dist_plot=ggplot(NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=10)) +
    labs(title=paste0("Pairwise Distances Between SNPs: ", 100*p, "% Tail"),x="Physical Distance (bp)", y = "Density") +
    geom_histogram(aes(x=as.vector(tail_dists), y = ..density.., col="Candidates", fill="Candidates"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists[[1]]), y = ..density.., col="Null 1", fill="Null 1"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists[[2]]), y = ..density.., col="Null 2", fill="Null 2"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists[[3]]), y = ..density.., col="Null 3", fill="Null 3"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists[[4]]), y = ..density.., col="Null 4", fill="Null 4"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists[[5]]), y = ..density.., col="Null 5", fill="Null 5"), alpha=0.1, bw=bw) +
    scale_discrete_manual(name='Data', values = c('Candidates' = 'red', 'Null 1' = 'blue', 'Null 2' = 'royalblue1', 'Null 3' = 'deepskyblue', 'Null 4' = 'darkcyan', 'Null 5' = 'cyan'), 
                          aesthetics = c("colour", "fill"))
  min_dist_plot=ggplot(NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=10)) +
    labs(title=paste0("Distance to Nearest SNP: ", 100*p, "% Tail"), x="Physical Distance (bp)", y = "Density") +
    geom_histogram(aes(x=as.vector(tail_dists_min), y = ..density.., col="Candidates", fill="Candidates"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists_min[[1]]), y = ..density.., col="Null 1", fill="Null 1"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists_min[[2]]), y = ..density.., col="Null 2", fill="Null 2"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists_min[[3]]), y = ..density.., col="Null 3", fill="Null 3"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists_min[[4]]), y = ..density.., col="Null 4", fill="Null 4"), alpha=0.1, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists_min[[5]]), y = ..density.., col="Null 5", fill="Null 5"), alpha=0.1, bw=bw) +
    scale_discrete_manual(name='Data', values = c('Candidates' = 'red', 'Null 1' = 'blue', 'Null 2' = 'royalblue1', 'Null 3' = 'deepskyblue', 'Null 4' = 'darkcyan', 'Null 5' = 'cyan'), 
                          aesthetics = c("colour", "fill"))
  for(plot in list(dist_plot, min_dist_plot)){print(plot)}
}

distance_between_cand_snps_2=function(xtx, p, n_samples=20, xmax=NA, bw=1){
  # Make sure there is a row per SNP
  xtx=unique(xtx[colnames(xtx)!='gene'])
  ## Data
  ### Select tail SNPs
  xtx.tail=xtx[xtx$empirical_p_med_cov<=p,]
  ## For each chr
  chrs=unique(xtx.tail$chr)
  tail_dists=c()
  tail_dists_min=c()
  null_dists=c()
  null_dists_min=c()
  for(chr in chrs[order(chrs)]){
    # Distance of tails SNPs
    tail_dists_chr=dist(xtx.tail[xtx.tail$chr==chr,'pos'])
    tail_dists=c(tail_dists, tail_dists_chr)
    ## Get minimum distances
    tail_dists_chr_m=as.matrix(tail_dists_chr)
    tail_dists_chr_m[tail_dists_chr_m==0]=NA
    tail_dists_chr_min=apply(tail_dists_chr_m, 1, FUN = min, na.rm = TRUE)
    tail_dists_min=c(tail_dists_min, tail_dists_chr_min)
    # Null SNPs
    for(i in 1:n_samples){
      ### Randomly sample the same number of SNPs in the tail from all sites
      xtx.null=xtx[sample(nrow(xtx), nrow(xtx.tail)),]
      # Distance of null SNPs
      null_dists_chr=dist(xtx.null[xtx.null$chr==chr,'pos'])
      null_dists=c(null_dists, null_dists_chr)
      ## Get minimum distances
      null_dists_chr_m=as.matrix(null_dists_chr)
      null_dists_chr_m[null_dists_chr_m==0]=NA
      null_dists_chr_min=apply(null_dists_chr_m, 1, FUN = min, na.rm = TRUE)
      null_dists_min=c(null_dists_min, null_dists_chr_min)
    }
  }
  # Plot
  dist_plot=ggplot(NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=10)) +
    xlim(0,xmax) +
    labs(title=paste0("Pairwise Distances Between SNPs: ", 100*p, "% Tail"),x="Physical Distance (bp)", y = "Density") +
    geom_histogram(aes(x=as.vector(tail_dists), y = ..density.., col="Candidates", fill="Candidates"), alpha=0.3, bw=bw) +
    geom_histogram(aes(x=as.vector(null_dists), y = ..density.., col="Null", fill="Null"), alpha=0.3, bw=bw) +
    scale_discrete_manual(name='Data', values = c('Candidates' = 'red', 'Null' = 'blue'), aesthetics = c("colour", "fill"))
  print(dist_plot)
}