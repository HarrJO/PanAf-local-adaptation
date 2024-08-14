# baypass_core_tools.R
## Tools for analysing the output from the BayPass core model
# Library
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plyr)
library(ape)
library(dendextend)
library(gplots)
library(cluster)
library(viridis) 
library(scatterpie)


plot.omega.PCA=function(omega, pops=NULL, title="PCA From Allele Frequency Covariance Matrix", subtitle=""){
  # "This function performs an eigen-decomposition of the scaled covariance matrix of the population allele frequencies (Ω in Figure 1) to allow representation in a two-dimension plot. This actually corresponds to a (between population) PCA–like analysis." 
  ## The above refers specifically to the function which can be accessed with `source("baypass_2.2/utils/baypass_utils.R")`. 
  ### I edit this function so I can make it pretty.
  if(is.null(pops)){pops=colnames(omega)}
  om.svd=svd(omega)
  eig=om.svd$d
  PCA=cbind(as.data.frame(pops), as.data.frame(om.svd$u))
  colnames(PCA)[1]="pop"
  # Calculate percentage variation for each PC
  pcent.var=100*eig/sum(eig)
  # Add subspecies column
  PCA$Subspecies <-""
  PCA[grepl('^c.', PCA$pop), 'Subspecies']="Central"
  PCA[grepl('^e.', PCA$pop), 'Subspecies']="Eastern"
  PCA[grepl('^n.', PCA$pop), 'Subspecies']="Nigeria-Cameroon"
  PCA[grepl('^w.', PCA$pop), 'Subspecies']="Western"
  # Group so I can plot it with the polygon connecting the points
  PC1.lab <- paste("PC1 (", round(pcent.var[1],digits = 2),"%)")
  PC2.lab <- paste("PC2 (", round(pcent.var[2],digits = 2),"%)")
  find_hull <- function(PCA) PCA[chull(PCA[,'V1'], PCA[,'V2']), ]
  hulls <- ddply(PCA, "Subspecies", find_hull)
  # Plot
  plot=ggplot(PCA, aes(x=V1, y=V2, fill=Subspecies, colour=Subspecies, shape=Subspecies, label = pop)) +
    theme_minimal()+ 
    geom_point() + 
    geom_polygon(data = hulls, alpha = 0.5) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5, size=20), plot.subtitle = element_text(hjust = 0.5, size=15),
          axis.title=element_text(size=15), axis.text = element_text(size=10),
          legend.title=element_text(size=15), legend.text=element_text(size=10)) +
    guides(fill=guide_legend(ncol=1)) + 
    labs(title=title, subtitle=subtitle, x=paste(PC1.lab), y=paste(PC2.lab)) + 
    geom_text_repel(size = 3, min.segment.length = unit(0.1, "lines"), show.legend = FALSE,max.overlaps = 15) +
    scale_fill_manual(values = c("Central"='green3', "Eastern"='darkorange', "Nigeria-Cameroon"='red', "Western"='blue')) +
    scale_color_manual(values = c("Central"='green3', "Eastern"='darkorange', "Nigeria-Cameroon"='red', "Western"='blue')) +
    scale_shape_manual(values = c("Central"=21, "Eastern"=24, "Nigeria-Cameroon"=22, "Western"=23)) + 
    coord_equal()
  #text(om.svd$u[,PC[1]],om.svd$u[,PC[2]],pop.names,col=col, cex=1, pos = 1)
  print(plot)
}


plot.omega.heatmap=function(omega, title="Allele Frequency\nCorrelation Matrix", out.file=NULL){
  # Plot the covraiance matrix
  par(oma=c(0,0,1.5,0))
  ### Set colour ramp
  #col=colorRampPalette(c("blue", "yellow", "red"))(100)
  col=plasma(100)
  ## As a heatmap and hierarchical clustering tree (using the average agglomeration method) 
  cor.mat=cov2cor(omega) 
  # Make distance matrix
  dist=as.dist(1 - cor.mat)
  # Make tree
  tree=hclust(dist, method="average")
  # Make dendrogram
  dend=as.dendrogram(tree)
  # Subspecies colours
  subsp.col=c()
  for(pop in rownames(cor.mat)){
    if(grepl('^c.', pop)==TRUE){subsp.col=c(subsp.col, "green3")}
    if(grepl('^e.', pop)==TRUE){subsp.col=c(subsp.col, "darkorange")}
    if(grepl('^n.', pop)==TRUE){subsp.col=c(subsp.col, "red")}
    if(grepl('^w.', pop)==TRUE){subsp.col=c(subsp.col, "blue")}
  }
  
  # Required to ensure column and row colours are in the correct order when using revC=T
  ## from https://stackoverflow.com/questions/34136096/r-gplots-heatmap-with-side-colours
  pdf(file = NULL) # This ensures it is not printed
  minimam.heatmap=heatmap.2(cor.mat,
                            Rowv = ladderize(dend),  
                            Colv = ladderize(dend), 
                            revC = FALSE,
                            dendrogram = "both", 
                            RowSideColors=subsp.col)
  dev.off()
  ordinary_order = minimam.heatmap$rowInd
  reversal = cbind(ordinary_order, rev(ordinary_order))
  rev_col = subsp.col[reversal[,2]]; rev_col = rev_col[order(reversal[,1])];
  
  if(!is.null(out.file)){pdf(file=out.file)}
  ## Actual plot
  cex=min(20/nrow(cor.mat), 1.5)
  heatmap.2(cor.mat,
            Rowv = ladderize(dend),  
            Colv = ladderize(dend), 
            revC = TRUE,
            dendrogram = "both", 
            col = col,
            trace = "none", density.info = "density",
            margins = c(15, 15),
            key.title=NA, key.xlab=NA, keysize=1,
            #key.ylab=NA, 
            main=title,
            ColSideColors=subsp.col, RowSideColors=rev_col,
            cexRow=cex, cexCol=cex)
  if(!is.null(out.file)){dev.off()}  
}

plot.omega.tree=function(omega, out.file=NULL){
  # Make correlation matrix
  cor.mat=cov2cor(omega)
  # Make dissimilarity matrix and convert to a phylo object
  tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
  #Plot
  plot.phylo(tree, 
             cex=0.5, font=2,
             main=expression("Hierarchical Clustering Tree Based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
  # Write to PDF
  if(!is.null(out.file)){
    pdf(file=out.file, pointsize=5, width=2.5, height=3.5)
    plot.phylo(tree, 
               cex=0.5, font=2,
               main=expression("Hierarchical Clustering Tree Based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
    dev.off()
  }
}

plot_xtx=function(xtx_exome, xtx_null, stat='XtXst.med', title="", alpha=0.2, binwidth=0.5, chi_df=NULL, xmin=0){
  if(stat=='XtXst.med'){stat_name="XtX*"}else{stat_name=stat}
  xtx_plot=ggplot(NULL)+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.title=element_text(size=15), axis.text = element_text(size=10),
          legend.title=element_text(size=20), legend.text=element_text(size=15), legend.position = c(0.9, 0.7), legend.key = element_rect(colour = "transparent")) +
    labs(title=paste0(title), x=stat_name, y="Density") +
    stat_bin(aes(x=xtx_exome[[stat]], y = ..density.., fill='Exome', colour='Exome'), alpha=1, binwidth=binwidth, geom="step",
             position=position_nudge(x=-0.5*binwidth)) +
    geom_histogram(aes(x=xtx_exome[[stat]], y = ..density.., fill='Exome', colour='Exome'), alpha=alpha, binwidth=binwidth, size=0)+
    xlim(xmin, NA)
  if(!is.null(chi_df)){
    xtx_plot=xtx_plot+
      stat_function(fun = dchisq, args = list(df = chi_df), aes(linetype = 'chi')) +
      scale_linetype_manual(name='', values = c('chi'="dashed"), labels = expression(chi^2))
  }
  if(!is.null(xtx_null)){
    xtx_plot=xtx_plot+
      stat_bin(aes(x=xtx_null[[stat]], y = ..density.., fill='Non-genic\n-chr21', colour='Non-genic\n-chr21'), alpha=1, binwidth=binwidth, geom="step",
               position=position_nudge(x=-0.5*binwidth)) +
      geom_histogram(aes(x=xtx_null[[stat]], y = ..density.., fill='Non-genic\n-chr21', colour='Non-genic\n-chr21'), alpha=alpha, binwidth=binwidth, size=0) +
      scale_discrete_manual(name='', values = c('Exome' = 'turquoise4', 'Non-genic\n-chr21' = 'violetred'), aesthetics = c("colour", "fill"))
  }else{
    xtx_plot=xtx_plot+
      scale_discrete_manual(name='', values = c('Exome' = 'turquoise4'), aesthetics = c("colour", "fill"))
  }
  print(xtx_plot)
}

# Assign FPRs
assign_fpr_core=function(xtx, null, n_bins=5, tail_bin_size=0.1){
  out=list()
  cover_bins=list()
  subsps=names(xtx)
  for(subsp in subsps){
    cat("-", subsp, "\n")
    # Assign SNPs to coverage bins
    ## Get all coverage values
    cover=xtx[[subsp]]$coverage
    ## Ensure that the lowest and highest bins contain 5% of the data points each
    low=max(head(cover[order(cover)], tail_bin_size*length(cover)))
    high=min(head(cover[order(-cover)], tail_bin_size*length(cover)))
    ## Separate into n bins between the low and high bin thresholds so we end up with n+2 bins over all (when you include the highest and lowest bins)
    cover_bins[[subsp]]=c(0, seq(from = low, to = high, by = (high-low)/(n_bins-2)), max(cover))
    ## Assign each SNP to a bin
    xtx[[subsp]]$coverage_bin=NA
    null[[subsp]]$coverage_bin=NA
    for(i in 1:(length(cover_bins[[subsp]])-1)){
      xtx[[subsp]][xtx[[subsp]]$coverage > cover_bins[[subsp]][i] & xtx[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
      null[[subsp]][null[[subsp]]$coverage > cover_bins[[subsp]][i] & null[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
    }
    # Loop over each bin and assign FDR values
    bins=unique(xtx[[subsp]]$coverage_bin)
    bin_fdr_out=data.frame(chr=numeric(), pos=numeric(), XtXst.med=numeric(), fdr=numeric(), empirical_p_med_cov=numeric())
    for(bin in bins[order(bins)]){
      cat("-- Bin:", bin, "/",length(bins),"\n")
      ## Extract SNPs from the bin in the null (and combine the information across random runs)
      ### Get null BFs 
      null_bin=null[[subsp]][(null[[subsp]]$coverage_bin==bin), c('chr', 'pos', 'XtXst.med')]
      null_bin$dataset='null'
      #### Get FPR values for each SNP
      null_bin$fpr=rank(-null_bin$XtXst.med, ties.method="max")/nrow(null_bin)
      #### Add real data
      xtx_bin=xtx[[subsp]][(xtx[[subsp]]$coverage_bin == bin), c('chr', 'pos', 'XtXst.med')]
      xtx_bin$dataset='real'
      xtx_bin$fpr=NA
      total_bin=rbind(null_bin, xtx_bin)
      #### Order dataset accordingly
      total_bin=total_bin[order(total_bin$XtXst.med, total_bin$dataset),]
      #### Loop over each row
      milestones=round(nrow(total_bin)*c(0.25, 0.5, 0.75))
      ##### The first value should have an FPR of 1
      total_bin[1, 'fpr']=1
      fpr_out=c(1)
      for(row in 2:nrow(total_bin)){
        if(row %in% milestones){
          cat("---", round(100*row/nrow(total_bin)),"% ")
        }
        if(is.na(total_bin[row, 'fpr'])){
          total_bin[row, 'fpr']=total_bin[row-1, 'fpr']
        }
      }
      cat("---", 100,"%\n")
      #total_bin$fpr=fpr_out
      #### Extract real data and add to file
      xtx_bin=total_bin[total_bin$dataset=='real',]
      ##### Calculate emprirical p values per bin too
      xtx_bin$empirical_p_med_cov=rank(-xtx_bin$XtXst.med, ties.method="max")/nrow(xtx_bin)
      bin_fdr_out=rbind(bin_fdr_out, xtx_bin)
    }
    # Add FDR values to dataframe    
    out[[subsp]]=merge(xtx[[subsp]][,c(c('chr', 'pos', 'XtXst.med', 'coverage', 'coverage_bin'))], bin_fdr_out, by=c('chr', 'pos', 'XtXst.med'), all.x=TRUE)
    out[[subsp]]=out[[subsp]][,colnames(out[[subsp]])!='dataset']
  }
  return(out)
}

plot_fpr_stats_core=function(xtx_fpr, null, fprs=c(0.005, 0.001, 0.0005), top.cov=0.025){
  subsps=names(xtx_fpr)
  fdr.df=NULL
  cat("Coverage Bin Stats\n")
  for(subsp in subsps){
    cat(subsp, "\n")
    for(fpr in c(1, fprs)){
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
                     'Threshold'= min(cands$`XtXst.med`),
                     'FDR'=fdr,
                     'True_Positives'=cand_n_snps*(1-fdr),
                     'False_Positives'=cand_n_snps*fdr)
      if(!is.null(fdr.df)){
        fdr.df=rbind(fdr.df, row)
      }else{
        fdr.df=row
      }
      # For each coverage bin
      bins=unique(xtx_fpr[[subsp]]$coverage_bin)
      for(bin in bins[order(bins)]){
        # Here I calculate the number of null SNPs which were used to calculate the FPR given the FPRs
        ## FPR is rank/total SNPs in bin - we want the maximum rank (which tells us how many null SNPs in that bin are over the threshold)
        ## First calculate total number of SNPs in bin by inputting rank 1 and and the minimum fpr (which must correspond to rank 1) 
        n_null_snps_bin=1/min(cands[cands$coverage_bin==bin,'fpr'])
        ## Now we have the total number of SNPs and the maximum FPR (corresponding to the SNP with the highest rank), we can calculate the highest rank
        ### Because in tied instances, rank takes the largest number, this tells us the total number of null SNPs above this threshold
        n_null_snps_bin_over_thresh=n_null_snps_bin*max(cands[cands$coverage_bin==bin,'fpr'])
        # Get the BF threshold for the bin
        thresh=min(cands[cands$coverage_bin==bin, 'XtXst.med'])
        # Add row for null data
        bin_df=rbind(bin_df, data.frame(fpr=fpr,
                                        bin=bin, 
                                        n_SNPs=n_null_snps_bin_over_thresh, 
                                        dataset='Null', 
                                        thresh=thresh))
        # Add row for real data (having it in long format makes it easier to plot)
        bin_df=rbind(bin_df, data.frame(fpr=fpr,
                                        bin=bin, 
                                        n_SNPs=nrow(cands[cands$coverage_bin==bin,]), 
                                        dataset='Exome', 
                                        thresh=thresh))
      }
      # Plot number of SNPs in bins over certain FPR
      n_snps_per_bin=ggplot(bin_df, aes(x=bin, y=n_SNPs, fill=dataset))+
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_bar(width=0.8, position = "dodge", stat="identity", col='black')+
        labs(title=paste0(subsp, ": Number of SNPs above FPR=",100*fpr,"% BF threshold"),
             x="Coverage bin", y="Count") +
        scale_discrete_manual(name="", 
                              labels = c("Null", "Exome"),
                              values = c('Null'='blue', 'Exome'='red'), 
                              aesthetics = "fill")
      print(n_snps_per_bin)
      # Plot BF thresholds per bin
      thresh_plot=ggplot(bin_df[bin_df$thresh!=Inf,], aes(x=bin, y=thresh))+
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_bar(width=0.8, position = "dodge", stat="identity")+
        labs(title=paste0(subsp, ": FPR=",100*fpr,"%"), y="XtXst Threshold", x="Coverage Bin") +
        coord_cartesian(ylim = c(min(19.5, 0.95*min(bin_df$thresh)), NA))
      print(thresh_plot)
      # Plot to check candidates do not have biased coverage
      ## Get (approximate) limits of bins using the maximum coverages observed in each bin
      cover_bins=aggregate(coverage ~ coverage_bin, xtx_fpr[[subsp]], function(x) max(x))$coverage
      ## Get xmax (to ignore the top x% of coverage values and make the graph readable)
      xmax=min(head(xtx_fpr[[subsp]][order(-xtx_fpr[[subsp]]$coverage),], nrow(xtx_fpr[[subsp]])*top.cov)$coverage)
      ## Parameters
      binwidth=100
      alpha=0.2
      lw=0.75
      ## Plot
      coverage_plot=ggplot(NULL) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
        #geom_density(data=xtx[[subsp]], aes(x=coverage, col="All Sites", fill="All Sites"), alpha=0.3, bw=100) +
        #geom_density(data=cands, aes(x=coverage, col="Candidates", fill="Candidates"), alpha=0.3, bw=100) +
        #geom_density(data=null[[subsp]], aes(x=coverage, col="Null", fill="Null"), alpha=0.3, bw=100) +
        stat_bin(data=xtx_fpr[[subsp]], aes(x=coverage, y=..density.., col="All Sites", fill="All Sites"), alpha=1, binwidth=binwidth, geom="step", size=lw,
                 position=position_nudge(x=-0.5*binwidth)) +
        geom_histogram(data=xtx_fpr[[subsp]], aes(x=coverage, y=..density.., fill="All Sites"), alpha=alpha, binwidth=binwidth, size=0) +
        stat_bin(data=cands, aes(x=coverage, y=..density.., col="Candidates", fill="Candidates"), alpha=1, binwidth=binwidth, geom="step", size=lw,
                 position=position_nudge(x=-0.5*binwidth)) +
        geom_histogram(data=cands, aes(x=coverage, y=..density.., fill="Candidates"), alpha=alpha, binwidth=binwidth, size=0) +
        stat_bin(data=null[[subsp]], aes(x=coverage, y=..density.., col="Null", fill="Null"), alpha=1, binwidth=binwidth, geom="step", size=lw,
                 position=position_nudge(x=-0.5*binwidth)) +
        geom_histogram(data=null[[subsp]], aes(x=coverage, y=..density.., fill="Null"), alpha=alpha, binwidth=binwidth, size=0) +
        geom_vline(xintercept = cover_bins, linetype='dashed', alpha=0.5) +
        labs(title=paste0(subsp, ", FPR=", fpr, ": Total Coverage Distribution"), 
             subtitle=paste0("Cropped X-axis to exclude the top ", 100*top.cov,"% coverages"), x="Total Coverage Across all Samples", y="Density") +
        scale_discrete_manual("", values = c('blue', 'red', 'green3'), aesthetics = c("colour", "fill")) + 
        xlim(0,xmax)
      print(coverage_plot)
    }
  }
  # Plot number of candidates and expected number
  cat("Number of Candidates\n")
  # For each fpr value
  for(fpr in fprs){
    # Plot number positives and null expectation for a specific FPR
    ## Select relevant columns 
    fdr.bar.df=fdr.df[fdr.df$FPR==fpr, c('Subspecies', 'Candidate_SNPs', 'Expected_SNPs')]
    if(nrow(fdr.bar.df)>0){
      ## Rename so plot looks nicer
      fdr.bar.df[fdr.bar.df$Subspecies=='all', 'Subspecies']="All Subspecies"
      fdr.bar.df[fdr.bar.df$Subspecies=='c', 'Subspecies']="Central"
      fdr.bar.df[fdr.bar.df$Subspecies=='e', 'Subspecies']="Eastern"
      fdr.bar.df[fdr.bar.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
      fdr.bar.df[fdr.bar.df$Subspecies=='n', 'Subspecies']="Nigeria-Cameroon"
      fdr.bar.df[fdr.bar.df$Subspecies=='w', 'Subspecies']="Western"
      #fdr.bar.df$Subspecies=factor(fdr.bar.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
      ## Plot
      ### Format file for input
      colnames(fdr.bar.df)=c('Subspecies', 'Number of Candidate SNPs', 'Expected Number of Candidate SNPs')
      fdr.bar.df$x='core'
      barplot=ggplot(fdr.bar.df, aes(x=x, y=`Number of Candidate SNPs`)) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              strip.text.x = element_text(face="bold"),
              strip.background = element_rect(fill="white"),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        geom_bar(stat="identity", fill='turquoise4', col='black') +
        #geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`, col="Null Expectation"),lty=2, size=1)+
        geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`), col="white", size=3)+
        geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`, col="Null Expectation"), size=2)+
        #geom_errorbar(aes(y = `Expected Number of Candidate SNPs`, ymin = `Expected Number of Candidate SNPs`, ymax = `Expected Number of Candidate SNPs`), col='red') + 
        facet_wrap(~Subspecies, ncol=length(unique(fdr.bar.df$Subspecies))) +
        labs(title=paste0("Number of Candidate SNPs"),
             subtitle=paste0("FPR = ", 100*fpr, "%"),
             x="") +
        scale_discrete_manual(name='', values = c("Null Expectation"='violetred'), aesthetics = c("col")) +
        scale_x_discrete(guide = guide_axis(angle = 45)) 
      print(barplot)
    }
  }
}

plot_fpr_stats_core.v2=function(xtx_fpr, null, fprs=c(0.005, 0.001, 0.0005), top.cov=0.025, N_SNPs=TRUE, BF_thresh=TRUE, coverage=TRUE){
  subsps=names(xtx_fpr)
  fdr.df=NULL
  cat("Coverage Bin Stats\n")
  for(subsp in subsps){
    if(subsp=='all'){
      subsp_name="All Subspecies"
    }else if(subsp=='ce'){
      subsp_name="Central-Eastern"
    }else if(subsp=='n'){
      subsp_name="Nigeria-Cameroon"
    }else if(subsp=='w'){
      subsp_name="Western"
    }else{
      subsp_name=subsp
    }
    for(fpr in c(fprs)){
      bin_df=data.frame(fpr=numeric(), bin=factor(), COVARIABLE_name=character(), n_SNPs=numeric(), dataset=character(), thresh=numeric())
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
                     'Threshold'= min(cands$`XtXst.med`),
                     'FDR'=fdr,
                     'True_Positives'=cand_n_snps*(1-fdr),
                     'False_Positives'=cand_n_snps*fdr)
      if(!is.null(fdr.df)){
        fdr.df=rbind(fdr.df, row)
      }else{
        fdr.df=row
      }
      # For each coverage bin
      bins=unique(xtx_fpr[[subsp]]$coverage_bin)
      for(bin in bins[order(bins)]){
        # Here I calculate the number of null SNPs which were used to calculate the FPR given the FPRs
        ## FPR is rank/total SNPs in bin - we want the maximum rank (which tells us how many null SNPs in that bin are over the threshold)
        ## First calculate total number of SNPs in bin by inputting rank 1 and and the minimum fpr (which must correspond to rank 1) 
        n_null_snps_bin=1/min(cands[cands$coverage_bin==bin,'fpr'])
        ## Now we have the total number of SNPs and the maximum FPR (corresponding to the SNP with the highest rank), we can calculate the highest rank
        ### Because in tied instances, rank takes the largest number, this tells us the total number of null SNPs above this threshold
        n_null_snps_bin_over_thresh=n_null_snps_bin*max(cands[cands$coverage_bin==bin,'fpr'])
        # Get the BF threshold for the bin
        thresh=min(cands[cands$coverage_bin==bin, 'XtXst.med'])
        # Add row for null data
        bin_df=rbind(bin_df, data.frame(fpr=fpr,
                                        bin=bin, 
                                        n_SNPs=n_null_snps_bin_over_thresh, 
                                        dataset='Null', 
                                        thresh=thresh))
        # Add row for real data (having it in long format makes it easier to plot)
        bin_df=rbind(bin_df, data.frame(fpr=fpr,
                                        bin=bin, 
                                        n_SNPs=nrow(cands[cands$coverage_bin==bin,]), 
                                        dataset='Exome', 
                                        thresh=thresh))
      }
      # Plot number of SNPs in bins over certain FPR
      if(N_SNPs){
        n_snps_per_bin=ggplot(bin_df, aes(x=bin, y=n_SNPs, fill=dataset))+
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_bar(width=0.8, position = "dodge", stat="identity", col='black')+
          labs(title=paste0(subsp_name, ", N SNPs at FPR<",100*fpr,"%"),
               x="Coverage bin", y="Count") +
          scale_discrete_manual(name="", 
                                #labels = c("Null", "Exome"),
                                values = c('Null'='green3', 'Exome'='red'), 
                                aesthetics = "fill")
        print(n_snps_per_bin)
      }
      
      # Plot BF thresholds per bin
      if(BF_thresh){
        thresh_plot=ggplot(bin_df[bin_df$thresh!=Inf,], aes(x=bin, y=thresh))+
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5),
                panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
          geom_bar(width=0.8, position = "dodge", stat="identity", col="black", fill="turquoise4")+
          labs(title=paste0(subsp_name, ", FPR<",100*fpr,"%"), y="XtX* threshold", x="Coverage bin") +
          coord_cartesian(ylim = c(0.95*min(bin_df$thresh), NA))
        print(thresh_plot)
      }
      # Plot to check candidates do not have biased coverage
      ## Get (approximate) limits of bins using the maximum coverages observed in each bin
      cover_bins=aggregate(coverage ~ coverage_bin, xtx_fpr[[subsp]], function(x) max(x))$coverage
      ## Get xmax (to ignore the top x% of coverage values and make the graph readable)
      xmax=min(head(xtx_fpr[[subsp]][order(-xtx_fpr[[subsp]]$coverage),], nrow(xtx_fpr[[subsp]])*top.cov)$coverage)
      ## Parameters
      binwidth=100
      alpha=0.2
      lw=0.75
      ## Plot
      if(coverage){
        coverage_plot=ggplot(NULL) + 
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
                legend.position = c(0.85, 0.85),
                legend.background = element_rect(fill='transparent')
                ) +
          #geom_density(data=xtx[[subsp]], aes(x=coverage, col="Exome", fill="Exome"), alpha=0.3, bw=100) +
          #geom_density(data=cands, aes(x=coverage, col="Candidates", fill="Candidates"), alpha=0.3, bw=100) +
          #geom_density(data=null[[subsp]], aes(x=coverage, col="Non-genic\n-chr21", fill="Non-genic\n-chr21"), alpha=0.3, bw=100) +
          stat_bin(data=xtx_fpr[[subsp]], aes(x=coverage, y=..density.., col="Exome", fill="Exome"), alpha=1, binwidth=binwidth, geom="step", size=lw,
                   position=position_nudge(x=-0.5*binwidth)) +
          geom_histogram(data=xtx_fpr[[subsp]], aes(x=coverage, y=..density.., fill="Exome"), alpha=alpha, binwidth=binwidth, size=0) +
          stat_bin(data=cands, aes(x=coverage, y=..density.., col="Candidates", fill="Candidates"), alpha=1, binwidth=binwidth, geom="step", size=lw,
                   position=position_nudge(x=-0.5*binwidth)) +
          geom_histogram(data=cands, aes(x=coverage, y=..density.., fill="Candidates"), alpha=alpha, binwidth=binwidth, size=0) +
          stat_bin(data=null[[subsp]], aes(x=coverage, y=..density.., col="Non-genic\n-chr21", fill="Non-genic\n-chr21"), alpha=1, binwidth=binwidth, geom="step", size=lw,
                   position=position_nudge(x=-0.5*binwidth)) +
          geom_histogram(data=null[[subsp]], aes(x=coverage, y=..density.., fill="Non-genic\n-chr21"), alpha=alpha, binwidth=binwidth, size=0) +
          geom_vline(xintercept = cover_bins, linetype='dashed', alpha=0.5) +
          labs(title=paste0(subsp_name, ", FPR<", 100*fpr, "%"), 
               subtitle=paste0("Cropped X-axis to exclude the top ", 100*top.cov,"% coverages"), x="Total coverage across samples", y="Density") +
          scale_discrete_manual("", values = c('darkorange', 'turquoise4', 'violetred'), aesthetics = c("colour", "fill")) + 
          xlim(0,xmax)
        print(coverage_plot)
      }
    }
  }
  # Plot number of candidates and expected number
  cat("Number of Candidates\n")
  # For each fpr value
  for(fpr in fprs){
    # Plot number positives and null expectation for a specific FPR
    ## Select relevant columns 
    fdr.bar.df=fdr.df[fdr.df$FPR==fpr, c('Subspecies', 'Candidate_SNPs', 'Expected_SNPs')]
    if(nrow(fdr.bar.df)>0){
      ## Rename so plot looks nicer
      fdr.bar.df[fdr.bar.df$Subspecies=='all', 'Subspecies']="All Subspecies"
      fdr.bar.df[fdr.bar.df$Subspecies=='c', 'Subspecies']="Central"
      fdr.bar.df[fdr.bar.df$Subspecies=='e', 'Subspecies']="Eastern"
      fdr.bar.df[fdr.bar.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
      fdr.bar.df[fdr.bar.df$Subspecies=='n', 'Subspecies']="Nigeria-Cameroon"
      fdr.bar.df[fdr.bar.df$Subspecies=='w', 'Subspecies']="Western"
      #fdr.bar.df$Subspecies=factor(fdr.bar.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
      ## Plot
      ### Format file for input
      colnames(fdr.bar.df)=c('Subspecies', 'Number of Candidate SNPs', 'Expected Number of Candidate SNPs')
      fdr.bar.df$x='core'
      barplot=ggplot(fdr.bar.df, aes(x=x, y=`Number of Candidate SNPs`)) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              strip.text.x = element_text(face="bold"),
              strip.background = element_rect(fill="white"),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        geom_bar(stat="identity", fill='turquoise4', col='black') +
        #geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`, col="Null Expectation"),lty=2, size=1)+
        geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`), col="white", size=3)+
        geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`, col="Null Expectation"), size=2)+
        #geom_errorbar(aes(y = `Expected Number of Candidate SNPs`, ymin = `Expected Number of Candidate SNPs`, ymax = `Expected Number of Candidate SNPs`), col='red') + 
        facet_wrap(~Subspecies, ncol=length(unique(fdr.bar.df$Subspecies))) +
        labs(title=paste0("Number of Candidate SNPs"),
             subtitle=paste0("FPR = ", 100*fpr, "%"),
             x="") +
        scale_discrete_manual(name='', values = c("Null Expectation"='violetred'), aesthetics = c("col")) +
        scale_x_discrete(guide = guide_axis(angle = 45)) 
      print(barplot)
    }
  }
}
plot_fpr_stats_core_main=function(xtx_fpr, null, fprs=c(0.005, 0.001, 0.0005), thresh_stat='`XtXst.med`'){
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
  bar_plot=ggplot(fdr.df, aes(x=FPR, y=Candidate_SNPs)) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.text.x = element_text(face="bold"),
          strip.background = element_rect(fill="white"),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    geom_bar(position="dodge", stat="identity", col='black', alpha=1, fill="turquoise4") +
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs), size = 1.4, col="white", width=0.85) + 
    geom_errorbar(aes(y = Expected_SNPs, ymin = Expected_SNPs, ymax = Expected_SNPs, col="Null\nexpectation"), size = 1, width=0.8) +
    facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies)), scales = "free") +
    #facet_wrap(~Subspecies, ncol=length(unique(fdr.df$Subspecies))) +
    labs(title=paste0("Number of candidate SNPs"),
         x="Esimated FPR threshold", y="Number of SNPs") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_discrete_manual(name='', values = c("Null\nexpectation"='violetred'), aesthetics = c("col"))
  print(bar_plot)
  return(fdr.df)
}

candidate_overlap_core=function(xtx, stat='empirical_p_med_cov', tails=c(0.005, 0.001, 0.0005), set_color=c("All Samples"="black", "Central-Eastern"="brown4", "Nigeria-Cameroon"="red", "Western"="blue")){
  subsps=names(xtx)[order(names(xtx))]
  set_color=set_color[set_color]
  if(stat=='empirical_p_med_cov'){
    stat_name="Empirical p"
  }else if(stat=='fpr'){
    stat_name="FPR"
  }else{
    stat_name=stat
    }
  # Overlap between subspecies datasets
  if(length(subsps>1)){
    for(tail in tails){
      tail_name=paste0(100*tail, "%")
      venn_list=list()
      for(subsp in subsps){
        #if(subsp %in% c('all', "All Samples", "All")){set_col=c(set_col, "black")}
        snps=xtx[[subsp]][xtx[[subsp]][[stat]]<=tail, c('chr', 'pos')]
        venn_list[[subsp]]=paste(snps$chr, snps$pos, sep="_")
      }
      if(length(venn_list)>1){
        venn=ggVennDiagram(venn_list) +
          theme(plot.title = element_text(hjust = 0.5, size=15), plot.subtitle = element_text(hjust = 0.5)) +
          labs(title=paste0("Candidate SNPs: ",stat_name," < ", tail_name)) +
          scale_fill_viridis() +
          #scale_color_manual(values = set_color) +
          scale_x_continuous(expand = expansion(mult = 0.2))
        print(venn)
        
      }
    }
    ## If the file has gene annotations, look at candidate gene overlap
    if('gene' %in% colnames(xtx[[1]])){
      for(tail in tails){
        tail_name=paste0(100*tail, "%")
        venn_list=list()
        for(subsp in subsps){
          genes=xtx[[subsp]][xtx[[subsp]][[stat]]<=tail, 'gene']
          venn_list[[subsp]]=genes
        }
        if(length(venn_list)>1){
          venn=ggVennDiagram(venn_list) +
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
            labs(title=paste0("Candidate Genes: ",stat_name," < ", tail_name)) +
            scale_fill_viridis() +
            #scale_color_manual(values = set_color) +
            scale_x_continuous(expand = expansion(mult = 0.2))
          print(venn)
        }
      }
    }
  }
}

plot_pi_prior=function(exome_a.b, null_a.b=NULL, subtitle=""){
  pi_prior_plot=ggplot() +
    theme_bw() +
    stat_function(fun = function(x) dbeta(x, exome_a.b[1], exome_a.b[2]), aes(color = "Exome")) +
    labs(title="Pi Prior Distribution", subtitle=subtitle, y="Density") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.text=element_text(size=10))
  if(is.null(null_a.b)){
    pi_prior_plot=pi_prior_plot+
      scale_discrete_manual(name='Data', values = c('Exome' = 'red'), aesthetics = c("colour"))
  }else{
    pi_prior_plot=pi_prior_plot+
      stat_function(fun = function(x) dbeta(x, null_a.b[1], null_a.b[2]), aes(color = "chr21")) +
      scale_discrete_manual(name='Data', values = c('Exome' = 'red', 'chr21' = 'blue'), aesthetics = c("colour"))
  }
  print(pi_prior_plot)
}

# Raw AF pi charts 
plot_raw_af_map=function(allele.freq_cand, chimp.ranges, title_prefix=""){
  # allele.freq_cand requires:
  ## a row per SNP per population
  ## Columns MRK, Population, M_P, Latitude and Longitude
  if(!all(c("MRK", "Population", "M_P", "Latitude", "Longitude") %in% colnames(allele.freq_cand))){
    stop("Dataframe must have columns: MRK, Population, M_P, Latitude and Longitude")}
  # chimp.ranges can be downloaded from https://www.iucnredlist.org/resources/spatial-data-download and chimp information extracted in R with:
  ## ranges=readOGR("../../maps/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
  ## chimp.ranges=ranges[ranges$binomial=="Pan troglodytes",]
  ## writeOGR(chimp.ranges, "../../maps/chimp_ranges/", "chimp_ranges", driver = "ESRI Shapefile")
  ## chimp.ranges=readOGR("../../../../maps/chimp_ranges/chimp_ranges.shp")
  
  # Make AAF column
  allele.freq_cand$DAF=allele.freq_cand$M_P
  allele.freq_cand$AAF=1-allele.freq_cand$DAF
  
  # New coordinates from plotting
  ## This is so they do not overlap eachother 
  ## They will be connected to their true location with a line
  allele.freq_cand$plot.long=allele.freq_cand$Longitude
  allele.freq_cand$plot.lat=allele.freq_cand$Latitude
  ### Adjust
  allele.freq_cand[allele.freq_cand$Population=='w.Kayan', 'plot.lat']=15
  allele.freq_cand[allele.freq_cand$Population=='w.Sangaredi', 'plot.lat']=10
  allele.freq_cand[allele.freq_cand$Population=='w.Sangaredi', 'plot.long']=-15
  allele.freq_cand[allele.freq_cand$Population=='w.Boe', 'plot.long']=-16
  allele.freq_cand[allele.freq_cand$Population=='w.Dindefelo', 'plot.long']=-13.5
  allele.freq_cand[allele.freq_cand$Population=='w.Dindefelo', 'plot.lat']=13.5
  allele.freq_cand[allele.freq_cand$Population=='w.Sobeya', 'plot.lat']=9.5
  allele.freq_cand[allele.freq_cand$Population=='w.Bafing', 'plot.lat']=14
  allele.freq_cand[allele.freq_cand$Population=='w.Bafing', 'plot.long']=-9
  allele.freq_cand[allele.freq_cand$Population=='w.BakounSobory', 'plot.long']=-10
  allele.freq_cand[allele.freq_cand$Population=='w.Grebo', 'plot.lat']=3.5
  allele.freq_cand[allele.freq_cand$Population=='w.Grebo', 'plot.long']=-8
  allele.freq_cand[allele.freq_cand$Population=='w.Djouroutou', 'plot.lat']=4
  allele.freq_cand[allele.freq_cand$Population=='w.Djouroutou', 'plot.long']=-6
  allele.freq_cand[allele.freq_cand$Population=='w.Tai', 'plot.long']=-6
  allele.freq_cand[allele.freq_cand$Population=='e.Nyungwe', 'plot.long']=28
  allele.freq_cand[allele.freq_cand$Population=='e.Nyungwe', 'plot.lat']=-3
  allele.freq_cand[allele.freq_cand$Population=='e.Gishwati', 'plot.long']=31
  allele.freq_cand[allele.freq_cand$Population=='e.Bwindi', 'plot.long']=28
  
  # Limits
  ## x
  xmin=min(allele.freq_cand$plot.long)
  xmax=max(allele.freq_cand$plot.long)
  xmin=xmin-(xmax-xmin)*0.1
  xmax=xmax+(xmax-xmin)*0.1
  ## x
  ymin=min(allele.freq_cand$plot.lat)
  ymax=max(allele.freq_cand$plot.lat)
  ymin=ymin-(ymax-ymin)*0.15
  ymax=ymax+(ymax-ymin)*0.15
  
  # Plot
  world_map=map_data("world")
  for(snp in unique(allele.freq_cand$MRK)){
    allele.freq_cand.snp=allele.freq_cand[allele.freq_cand$MRK==snp,]
    ## Title
    if('chr' %in% colnames(allele.freq_cand.snp) & 'pos' %in% colnames(allele.freq_cand.snp)){
      title=paste0(title_prefix, "Raw Allele Frequencies, chr", unique(allele.freq_cand.snp$chr),":", unique(allele.freq_cand.snp$pos))
    }else{
      title=paste0(title_prefix, "Raw Allele Frequencies, SNP: ", snp)
    }
    plot=ggplot(allele.freq_cand.snp, aes(x=plot.long, y=plot.lat)) +
      # Map
      theme_void() +
      geom_polygon(data=world_map, aes(x = long, y = lat, group = group), fill="grey93", colour = "darkgrey") +
      coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) +
      # Subspecies ranges
      geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="troglodytes",], aes(x = long, y = lat, group = group), fill="green3", alpha=0.3) +
      geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="schweinfurthii",], aes(x = long, y = lat, group = group), fill="darkorange", alpha=0.3) +
      geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="ellioti",], aes(x = long, y = lat, group = group), fill="red", alpha=0.3) +
      geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="verus",], aes(x = long, y = lat, group = group), fill="blue", alpha=0.3) +
      # Pie charts
      geom_segment(aes(x = Longitude, y = Latitude, xend = plot.long, yend = plot.lat)) +
      geom_text_repel(aes(label=Population), 
                      size = 2, show.legend = FALSE, box.padding=0.25, point.padding=5, nudge_y=1, segment.alpha=0) +
      # General
      labs(title=title, x="", y = "") +
      theme(plot.title = element_text(hjust = 0.5))
    print(plot)
  }
}

# Standradised allele frequencies 
plot_af_map=function(allele.freq_cand, af_stat="M_P", title_prefix="", fill_col=NULL){
  # allele.freq_cand requires:
  ## a row per SNP per population
  ## Columns MRK, Population, M_Pstd, Latitude and Longitude
  if(!all(c("MRK", "Population", "M_Pstd", "Latitude", "Longitude") %in% colnames(allele.freq_cand))){
    stop("Dataframe must have columns: MRK, Population, M_Pstd, Latitude and Longitude")}
  # Load maps
  if(!exists("world")){
    world_map=map_data("world")
  }
  if(!exists("chimp.ranges")){
    chimp.ranges=readOGR("~/OneDrive - University College London/Projects/PanAf/maps/chimp_ranges/chimp_ranges.shp")
  }
  ## Downloaded version 5.0.0 from https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-rivers-lake-centerlines/
  if(!exists("rivers")){
    rivers=readOGR("~/OneDrive - University College London/Projects/PanAf/maps/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp")
    rivers=rivers[grepl("Congo|Niger|Ubangi|Ogooué|Tanganyika|Sanaga|Benue|Uelé|Lualaba|Chambeshi|Chinko|Mbomou|Uele", rivers$name_en),]
  }
  ## Download water bodies Zip (last updated Nov 14, 2018) https://datacatalog.worldbank.org/search/dataset/0040797
  if(!exists("lakes")){
    lakes=readOGR("~/OneDrive - University College London/Projects/PanAf/maps/africawaterbody/Africa_waterbody.shp")
  }
  # New coordinates from plotting
  ## This is so they do not overlap eachother 
  ## They will be connected to their true location with a line
  allele.freq_cand$plot.long=allele.freq_cand$Longitude
  allele.freq_cand$plot.lat=allele.freq_cand$Latitude
  ### Adjust
  allele.freq_cand[allele.freq_cand$Population=='w.Kayan', 'plot.lat']=15
  allele.freq_cand[allele.freq_cand$Population=='w.Sangaredi', 'plot.lat']=10
  allele.freq_cand[allele.freq_cand$Population=='w.Sangaredi', 'plot.long']=-15
  allele.freq_cand[allele.freq_cand$Population=='w.Boe', 'plot.long']=-16
  allele.freq_cand[allele.freq_cand$Population=='w.Dindefelo', 'plot.long']=-13.5
  allele.freq_cand[allele.freq_cand$Population=='w.Dindefelo', 'plot.lat']=13.5
  allele.freq_cand[allele.freq_cand$Population=='w.Sobeya', 'plot.lat']=9.5
  allele.freq_cand[allele.freq_cand$Population=='w.Bafing', 'plot.lat']=14
  allele.freq_cand[allele.freq_cand$Population=='w.Bafing', 'plot.long']=-9
  allele.freq_cand[allele.freq_cand$Population=='w.BakounSobory', 'plot.long']=-10
  allele.freq_cand[allele.freq_cand$Population=='w.Grebo', 'plot.lat']=3.5
  allele.freq_cand[allele.freq_cand$Population=='w.Grebo', 'plot.long']=-8
  allele.freq_cand[allele.freq_cand$Population=='w.Djouroutou', 'plot.lat']=4
  allele.freq_cand[allele.freq_cand$Population=='w.Djouroutou', 'plot.long']=-6
  allele.freq_cand[allele.freq_cand$Population=='w.Tai', 'plot.long']=-6
  allele.freq_cand[allele.freq_cand$Population=='e.Nyungwe', 'plot.long']=28
  allele.freq_cand[allele.freq_cand$Population=='e.Nyungwe', 'plot.lat']=-3
  allele.freq_cand[allele.freq_cand$Population=='e.Gishwati', 'plot.long']=31
  allele.freq_cand[allele.freq_cand$Population=='e.Bwindi', 'plot.long']=28
  # If were doing pie charts
  if(af_stat=="M_P"){
    # Make AAF column
    allele.freq_cand$DAF=allele.freq_cand$M_P
    allele.freq_cand$AAF=1-allele.freq_cand$DAF
  }
  # Limits
  ## x
  xmin=min(allele.freq_cand$plot.long)
  xmax=max(allele.freq_cand$plot.long)
  xmin=xmin-(xmax-xmin)*0.1
  xmax=xmax+(xmax-xmin)*0.1
  ## x
  ymin=min(allele.freq_cand$plot.lat)
  ymax=max(allele.freq_cand$plot.lat)
  ymin=ymin-(ymax-ymin)*0.15
  ymax=ymax+(ymax-ymin)*0.15
  
  # Point size
  point.size=100/(ymax-ymin)
  # Plot
  world_map=map_data("world")
  
  for(snp in unique(allele.freq_cand$chr_pos)){
    allele.freq_cand.snp=allele.freq_cand[allele.freq_cand$chr_pos==snp,]
    ## Title
    if('chr' %in% colnames(allele.freq_cand.snp) & 'pos' %in% colnames(allele.freq_cand.snp)){
      title=paste0(title_prefix, ", chr", unique(allele.freq_cand.snp$chr),":", unique(allele.freq_cand.snp$pos))
    }else{
      title=paste0(title_prefix, ", SNP: ", snp)
    }
    plot=ggplot(allele.freq_cand.snp, aes(x=plot.long, y=plot.lat)) +
      ## Map
      theme_void() +
      #theme(panel.background = element_rect(fill = 'cadetblue1')) +
      geom_polygon(data=world_map, aes(x = long, y = lat, group = group), fill="navajowhite3", colour = "navajowhite4", alpha=0.25) +
      coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) +
      labs(title=title) +
      theme(plot.title = element_text(hjust = 0.5, size=15)) +
      ## Subspecies ranges
      geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="troglodytes",], aes(x = long, y = lat, group = group), fill="green3", alpha=0.4) +
      geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="schweinfurthii",], aes(x = long, y = lat, group = group), fill="darkorange", alpha=0.4) +
      geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="ellioti",], aes(x = long, y = lat, group = group), fill="red", alpha=0.4) +
      geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="verus",], aes(x = long, y = lat, group = group), fill="blue", alpha=0.4) +
      ## Water
      geom_polygon(data=lakes, aes(x = long, y = lat, group=group), fill="deepskyblue", alpha=0.5)+
      geom_path(data=rivers, aes(x = long, y = lat, group = group), col="deepskyblue", alpha=0.5) +
      # Lines connecting points and names
      geom_segment(aes(x = Longitude, y = Latitude, xend = plot.long, yend = plot.lat)) +
      # General
      labs(title=title, x="", y = "") +
      theme(plot.title = element_text(hjust = 0.5))
    if(af_stat=="M_Pstd"){
      plot=plot+
        geom_point(data=allele.freq_cand.snp, aes(x=plot.long, y=plot.lat, fill=M_Pstd), size=point.size, pch=21, colour='black') +
        scale_fill_gradient2(low = "blue", mid = "white", midpoint = 0, high = "red", space = "Lab" )
    }else if(af_stat=="M_P"){
      plot=plot+
        geom_scatterpie(data=allele.freq_cand.snp, aes(x=plot.long, y=plot.lat, group=Population, r=0.75), 
                        cols=c("DAF","AAF"), legend_name = "Alleles") +
        scale_discrete_manual(name='Alleles', values = c('DAF' = 'black', 'AAF' = 'white'), aesthetics = c("fill"))
    }else{
      stop("'af_stat' must be M_P or M_Pstd")
    }
    # Plot names
    if(!is.null(fill_col)){
      plot=plot+
        geom_text_repel(aes_string(label='Population', col=fill_col), 
                        #col='black',
                        fontface= "bold",
                        min.segment.length = unit(9999, 'lines'),
                        size=3,
                        seed=2,
                        bg.color = "white",
                        bg.r = 0.1,
                        max.overlaps = Inf) +
        #scale_colour_gradientn(name=fill_col, colours = colorRampPalette(c('goldenrod1', 'orange3', 'yellow4', 'darkgreen'))(100), aesthetics = c("colour"))
        scale_colour_gradientn(name=fill_col, colours = colorRampPalette(c('goldenrod1', 'khaki4', 'darkgreen'))(100), aesthetics = c("colour"))
    }else{
      plot=plot+
        geom_text_repel(aes(label=Population), 
                        col='black',
                        fontface= "bold",
                        min.segment.length = unit(9999, 'lines'),
                        size=3,
                        seed=2,
                        bg.color = "white",
                        bg.r = 0.1,
                        max.overlaps = Inf)
    }
    print(plot)
  }
}




