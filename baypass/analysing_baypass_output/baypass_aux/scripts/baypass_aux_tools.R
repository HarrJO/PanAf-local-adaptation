# baypass_aux_tools.R
## Tools for analysing the output from the BayPass AUX model
# Library
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(ggVennDiagram)
library(grid)
library(ggrepel)
library(viridis)
library(rgdal)
library(wesanderson)

# Plot
# Volcano Plot
aux_volcano=function(betai, cov_cols=NULL){
  if(is.null(cov_cols)){
    # This function just selects equally spaced hues around the colour wheel (which is ggplot's default)
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cov_cols = gg_color_hue(length(unique(betai$COVARIABLE_name)))
  }
  # c(t, r, b, l)
  # Beta density
  ## pseudo_log_trans() is great, removes the issue of 0s completely messing up the axis scale
  x_plot=ggplot(betai, aes_string(x='M_Beta', col='COVARIABLE_name')) +
    theme_bw() +
    geom_density(alpha=0.1) +
    scale_discrete_manual(name='Covariable', values = cov_cols, aesthetics = c("colour", "fill")) +
    theme(axis.title.x=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "mm"),
          legend.justification = "top",
          legend.key.size = unit(min(23/length(unique(betai$COVARIABLE_name)), 6), 'mm'),
          legend.title = element_text(size=min(45/length(unique(betai$COVARIABLE_name)), 15)),
          legend.text = element_text(size=min(40/length(unique(betai$COVARIABLE_name)), 12))) +
    labs(y=expression(log[10]~(Density))) + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), expand = c(0, 0))
  ## Extract legend
  legend=get_legend(x_plot)
  ## Remove legend from plot
  x_plot=x_plot + theme(legend.position = "none")
  
  # Beta vs BF
  main_plot=ggplot(data=betai, aes(x=M_Beta, y=`BF(dB)`, col=COVARIABLE_name, fill=COVARIABLE_name)) + 
    theme_bw() +
    geom_point(alpha=0.05, shape=1) +
    labs(title=paste0("")) +
    scale_discrete_manual(name='Covariable', values = cov_cols, aesthetics = c("colour", "fill")) +
    theme(legend.position = "none", 
          plot.title=element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "mm"))
  
  # BF density
  y_plot=ggplot(betai, aes_string(x='`BF(dB)`', col='COVARIABLE_name')) +
    theme_bw() +
    geom_density(alpha=0.1) +
    scale_discrete_manual(name='Covariable', values = cov_cols, aesthetics = c("colour", "fill")) +
    theme(axis.title.y=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0, 0, 0), "mm")) +
    labs(y="Density") +
    coord_flip()
  
  # Plot
  ## The following is required to ensure axes are lined up
  ## No need to perform proper alignment (infact it mess things up as you try to align with the legend axis) because we are plotting the exact same data on the axes so simply matching the lengths will line them up perfectly.
  plots <- list(x_plot, legend, main_plot, y_plot)
  grobs <- lapply(plots, as_grob)
  grobs[[1]]$widths=grobs[[3]]$widths
  grobs[[4]]$heights=grobs[[3]]$heights
  vol_plot=plot_grid(plotlist = grobs, ncol = 2, rel_widths = c(3, 1), rel_heights = c(1, 3))
  print(vol_plot)
}

aux_volcano.v2=function(betai, x='M_Beta.median', y='`BF(dB).median`', colour_col='COVARIABLE_name', colours=NULL, log_density=TRUE, alpha=NULL){
  if(x=='M_Beta.median'){x_lab=expression(Correlation~coefficent~(beta))}else{x_lab=x}
  if(y=='`BF(dB).median`'){y_lab="Bayes factor (dB)"}else{y_lab=y}
  if(is.null(colours)){
    # This function just selects equally spaced hues around the colour wheel (which is ggplot's default)
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    colours = gg_color_hue(length(unique(betai[[colour_col]])))
  }
  # c(t, r, b, l)
  # Beta density
  ## pseudo_log_trans() is great, removes the issue of 0s completely messing up the axis scale
  n_bins=51
  binwidth_x=(max(betai[[x]])-min(betai[[x]]))/n_bins
  binwidth_y=(max(betai[[y]])-min(betai[[y]]))/n_bins
  x_plot=ggplot(betai, aes_string(x=x, col=colour_col, y='..density..')) +
    theme_bw() +
    #geom_density(alpha=0.1) +
    stat_bin(bins=n_bins, geom="step", position=position_nudge(x=-0.5*binwidth_x)) +
    scale_discrete_manual(name='', values = colours, aesthetics = c("colour", "fill")) +
    theme(axis.title.x=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "mm"),
          legend.justification = "top",
          legend.key.size = unit(min(23/length(unique(betai[[colour_col]])), 6), 'mm'),
          legend.title = element_text(size=min(45/length(unique(betai[[colour_col]])), 15)),
          legend.text = element_text(size=min(40/length(unique(betai[[colour_col]])), 10))) +
    guides(colour = guide_legend(override.aes = list(alpha=1))) +
    xlim(min(betai[[x]]), max(betai[[x]]))
  ## Extract legend
  legend=get_legend(x_plot)
  ## Remove legend from plot
  x_plot=x_plot + theme(legend.position = "none")
  
  # Beta vs BF
  ## If the file is huge I take a random subset to plot the points to speed things up and also no points appear when trying to plot loads
  betai_points=betai[sample(nrow(betai), min(10000000, nrow(betai))),]
  if(is.null(alpha)){
    alpha=max(0.01, min(c(0.5, 40000/nrow(betai_points))))
  }
  
  main_plot=ggplot(data=betai_points, aes_string(x=x, y=y, col=colour_col, fill=colour_col)) + 
    theme_bw() +
    geom_point(alpha=alpha, shape=4) +
    labs(title=paste0(""), x=x_lab, y=y_lab) +
    scale_discrete_manual(name='', values = colours, aesthetics = c("colour", "fill")) +
    theme(plot.title=element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "mm")) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size = 5, shape=15)))
  ## Extract large legend (in case we need to plot separately)
  big_legend=get_legend(main_plot)
  main_plot=main_plot+theme(legend.position = "none")
  # BF density
  y_plot=ggplot(betai, aes_string(x=y, col=colour_col, y='..density..')) +
    theme_bw() +
    #geom_density(alpha=0.1) +
    stat_bin(bins=n_bins, geom="step", position="identity") +
    scale_discrete_manual(name='', values = colours, aesthetics = c("colour", "fill")) +
    theme(axis.title.y=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0, 0, 0), "mm")) +
    coord_flip()
  
  if(log_density){
    #x_plot=x_plot + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), expand = c(0, 0)) + labs(y=expression(log[10]~(Density)))
    x_plot=x_plot + scale_y_log10() + labs(y=expression(log[10]~density)) +
      theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
    #y_plot=y_plot + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), expand = c(0, 0)) + labs(y=expression(log[10]~(Density)))
    y_plot=y_plot + scale_y_log10() + labs(y=expression(log[10]~density)) +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  }else{
    x_plot=x_plot + labs(y="Density")
    y_plot=y_plot + labs(y="Density")
  }
  
  # Plot
  ## The following is required to ensure axes are lined up
  ## No need to perform proper alignment (in fact it mess things up as you try to align with the legend axis) because we are plotting the exact same data on the axes so simply matching the lengths will line them up perfectly.
  ## If there are more than 10 categories we will print the legend as a seperate plot
  n_cat=length(unique(betai[[colour_col]]))
  if(n_cat>10){
    plots=list(x_plot, NULL, main_plot, y_plot)
  }else{
    plots=list(x_plot, legend, main_plot, y_plot)
  }
  grobs <- lapply(plots, as_grob)
  grobs[[1]]$widths=grobs[[3]]$widths
  grobs[[4]]$heights=grobs[[3]]$heights
  vol_plot=plot_grid(plotlist = grobs, ncol = 2, rel_widths = c(3, 1), rel_heights = c(1, 3))
  print(vol_plot)
  if(n_cat>10){
    grid.newpage()
    grid.draw(big_legend)
  }
}

# Plot allele frequencies of candidates in different habitat types
env_vs_AF=function(betai, allele.freq, pop_cov, stat="M_Pstd", direction='both', bf.thresh=20, beta.thresh=0.1){
  cov=unique(betai$COVARIABLE_name)
  # Check there is only results for one covariable in the Betai dataframe
  if(length(cov)>1){stop('Supply Betai dataframe with results from only one covariable')}
  # Select candidates (depending on whether beta is +ve or -ve)
  if(direction=='positive'){
    cand_betai=betai[(betai$`BF(dB)`>bf.thresh) & (betai$M_Beta>beta.thresh),]
    subtitle=paste0("BF > ",bf.thresh, " & Beta > ", beta.thresh)
  }
  if(direction=='negative'){
    cand_betai=betai[(betai$`BF(dB)`>bf.thresh) & (betai$M_Beta<-beta.thresh),]
    subtitle=paste0("BF > ",bf.thresh, " & Beta < -", beta.thresh)
  }
  if(direction=='both'){
    cand_betai=betai[(betai$`BF(dB)`>bf.thresh) & (abs(betai$M_Beta)>beta.thresh),]
    subtitle=paste0("BF > ",bf.thresh, " & abs(Beta) > ", beta.thresh)
  }
  if(!(direction %in% c('positive', 'negative', 'both'))){stop("Select beta candidates (either 'positive', 'negative' or 'both')")}
  # If there are any candidates...
  n_cand=nrow(cand_betai)
  if(n_cand>0){
    # Select allele frequencies for candidate sites
    cand_allele.freq=allele.freq[allele.freq$MRK %in% cand_betai$MRK,]
    # Add covariable information
    cand_allele.freq=merge(cand_allele.freq, pop_cov, by='Population')
    # Add Beta values
    cand_allele.freq=merge(cand_allele.freq, unique(betai[, c('MRK', 'M_Beta')]), all.x=TRUE)
    # Get Beta direction (to colour lines and boxplots)
    cand_allele.freq$Direction=NA
    cand_allele.freq[cand_allele.freq$M_Beta<0, 'Direction']='Negative'
    cand_allele.freq[cand_allele.freq$M_Beta>0, 'Direction']='Positive'
    # Get the average allele frequency for a given SNP in a given habitat (and keep all the relevant columns which are specific to a particular SNP anyway)
    cand_allele.freq=aggregate(cand_allele.freq[[stat]], 
                               by = list(cand_allele.freq$MRK, cand_allele.freq$Vrishti_habitats, cand_allele.freq$M_Beta, cand_allele.freq$Direction), FUN = mean)
    colnames(cand_allele.freq)=c('MRK', 'Vrishti_habitats', 'M_Beta', 'Direction', stat)
    # Plot
    cand_plot=ggplot(cand_allele.freq, aes_string(x='Vrishti_habitats', y=stat, col='Direction', fill='factor(Direction)')) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      labs(title=paste0(cov, " Candidates"), 
           subtitle=paste0(subtitle, " (", n_cand, " SNPs)"),
           y=paste0(stat, " (mean across populations)")) +
      geom_line(aes(group=MRK), alpha=min(0.25, 150/n_cand)) +
      geom_boxplot(width=0.3, outlier.size=0.5, outlier.shape=3) +
      scale_discrete_manual(name='Direction', values = c('Negative' = 'red', 'Positive'='green3'), aesthetics = c("colour")) +
      scale_discrete_manual(name='Direction', values = c('Negative' = 'darkred', 'Positive'='darkgreen'), aesthetics = c("fill"))
    print(cand_plot)
  }else(cat("No candidates for", unique(betai$COVARIABLE_name), ": ", subtitle, "\n"))
}

env_vs_AF.v2=function(betai, allele.freq, pop_cov, stat="M_Pstd", direction='both', subtitle='', beta='M_Beta'){
  cov=unique(betai$COVARIABLE_name)
  # Add chr_pos column
  betai$chr_pos=paste(betai$chr, betai$pos, sep="_")
  allele.freq$chr_pos=paste(allele.freq$chr, allele.freq$pos, sep="_")
  # Get column name
  env_column=colnames(pop_cov)[colnames(pop_cov)!='Population']
  # Check there is only results for one covariable in the Betai dataframe
  if(length(cov)>1){stop('Supply Betai dataframe with results from only one covariable')}
  # Select candidates (depending on whether beta is +ve or -ve)
  if(direction=='positive'){
    betai=betai[betai[[beta]]>0,]
  }
  if(direction=='negative'){
    betai=betai[betai[[beta]]<0,]
  }
  if(!(direction %in% c('positive', 'negative', 'both'))){stop("Select beta candidates (either 'positive', 'negative' or 'both')")}
  # If there are any candidates...
  n_cand=nrow(betai)
  if(n_cand>0){
    # Select allele frequencies for candidate sites
    allele.freq.select=allele.freq[allele.freq$chr_pos %in% betai$chr_pos,]
    # Add covariable information
    allele.freq.select=merge(allele.freq.select, pop_cov, by='Population')
    # Add Beta values
    allele.freq.select=merge(allele.freq.select, unique(betai[, c('chr_pos', beta)]), all.x=TRUE)
    # Get Beta direction (to colour lines and boxplots)
    allele.freq.select$Direction=NA
    allele.freq.select[allele.freq.select[[beta]]<0, 'Direction']='Negative'
    allele.freq.select[allele.freq.select[[beta]]>0, 'Direction']='Positive'
    # Get the average allele frequency for a given SNP in a given habitat (and keep all the relevant columns which are specific to a particular SNP anyway)
    allele.freq.select=aggregate(allele.freq.select[[stat]], 
                               by = list(allele.freq.select$chr_pos, allele.freq.select[[env_column]], allele.freq.select[[beta]], allele.freq.select$Direction), FUN = mean)
    colnames(allele.freq.select)=c('chr_pos', env_column, beta, 'Direction', stat)
    # Plot
    cand_plot=ggplot(allele.freq.select, aes_string(x=env_column, y=stat, col='Direction', fill='factor(Direction)')) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=25), plot.subtitle = element_text(hjust = 0.5, size=15)) +
      labs(title=paste0(cov, " Candidates"), 
           subtitle=paste0(subtitle, " (", n_cand, " SNPs)"),
           y=paste0(stat, " (mean across populations)")) +
      geom_line(aes(group=chr_pos), alpha=min(0.25, 150/n_cand)) +
      geom_boxplot(width=0.3, outlier.size=0.5, outlier.shape=3) +
      scale_discrete_manual(name='Direction', values = c('Negative' = 'red', 'Positive'='green3'), aesthetics = c("colour")) +
      scale_discrete_manual(name='Direction', values = c('Negative' = 'darkred', 'Positive'='darkgreen'), aesthetics = c("fill"))
    ## If dealing with raw allele frequencies, set y limits to 0,1
    if(stat=="M_P"){cand_plot=cand_plot + ylim(0,1)}
    ## Print plot
    print(cand_plot)
  }else(cat("No sites with ",direction," beta for", unique(betai$COVARIABLE_name), "\n"))
}

env_vs_AF.v3=function(betai, allele.freq, pop_cov, stat="M_Pstd", direction='both', subtitle='', beta='M_Beta.median'){
  cov=unique(betai$COVARIABLE_name)
  # Add chr_pos column
  betai$chr_pos=paste(betai$chr, betai$pos, sep="_")
  allele.freq$chr_pos=paste(allele.freq$chr, allele.freq$pos, sep="_")
  # Get column name
  env_column=colnames(pop_cov)[colnames(pop_cov)!='Population']
  # Check there is only results for one covariable in the Betai dataframe
  if(length(cov)>1){stop('Supply Betai dataframe with results from only one covariable')}
  # Select candidates (depending on whether beta is +ve or -ve)
  if(direction=='positive'){
    betai=betai[betai[[beta]]>0,]
  }
  if(direction=='negative'){
    betai=betai[betai[[beta]]<0,]
  }
  if(!(direction %in% c('positive', 'negative', 'both'))){stop("Select beta candidates (either 'positive', 'negative' or 'both')")}
  # If there are any candidates...
  n_cand=nrow(betai)
  if(n_cand>0){
    # Select allele frequencies for candidate sites
    allele.freq.select=allele.freq[allele.freq$chr_pos %in% betai$chr_pos,]
    # Add covariable information
    allele.freq.select=merge(allele.freq.select, pop_cov, by='Population')
    # Add Beta values
    allele.freq.select=merge(allele.freq.select, unique(betai[, c('chr_pos', beta)]), all.x=TRUE)
    # Get Beta direction (to colour lines and boxplots)
    allele.freq.select$Direction=NA
    if(direction!='negative'){allele.freq.select[allele.freq.select[[beta]]<0, 'Direction']='Negative'}
    if(direction!='positive'){allele.freq.select[allele.freq.select[[beta]]>0, 'Direction']='Positive'}
    # Get the average allele frequency for a given SNP in a given habitat (and keep all the relevant columns which are specific to a particular SNP anyway)
    allele.freq.select=aggregate(allele.freq.select[[stat]], 
                                 by = list(allele.freq.select$chr_pos, allele.freq.select[[env_column]], allele.freq.select[[beta]], allele.freq.select$Direction), FUN = mean)
    colnames(allele.freq.select)=c('chr_pos', env_column, beta, 'Direction', stat)
    # Plot
    cand_plot=ggplot(allele.freq.select, aes_string(x=env_column, y=stat, col='Direction', fill='factor(Direction)')) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=25), plot.subtitle = element_text(hjust = 0.5, size=15)) +
      labs(title=paste0(cov, " :", stat), 
           subtitle=paste0(subtitle, " (", n_cand, " SNPs)"),
           y=paste0(stat, " (mean across populations)")) +
      geom_line(aes(group=chr_pos), alpha=min(0.25, 150/n_cand)) +
      scale_discrete_manual(name='Direction', values = c('Negative' = 'red', 'Positive'='green3'), aesthetics = c("colour")) +
      scale_discrete_manual(name='Direction', values = c('Negative' = 'darkred', 'Positive'='darkgreen'), aesthetics = c("fill"))
    ## If we are dealin with discrete habitats then add boxplots
    if(env_column %in% c('aleman_hab', 'aleman_hab_1cat')){
      cand_plot=cand_plot+geom_boxplot(width=0.3, outlier.size=0.5, outlier.shape=3)
    }
    ## If dealing with raw allele frequencies, set y limits to 0,1
    if(stat=="M_P"){cand_plot=cand_plot + ylim(0,1)}
    ## Print plot
    print(cand_plot)
  }else(cat("No sites with ",direction," beta for", unique(betai$COVARIABLE_name), "\n"))
}

env_vs_AF.v4=function(betai, allele.freq, pop_cov, env_column, stat="M_Pstd", direction='both', title='', 
                      beta='M_Beta.median', pop_colour_column=NULL, cor_method='pearson', big_line=NULL, 
                      colour_column='Direction', legend=TRUE, pop_labels=TRUE){
  if(!colour_column %in% c('Direction', colnames(betai))){
    stop("colour_column should be either 'Direction' (default) or another column in the betai file")
  }
  if(stat=='M_P'){
    stat_name='Unstandardised allele frequency'
  }else if(stat=='M_Pstd'){
     stat_name='Standardised allele frequency'
  }else if(stat=='raw_AF'){
    stat_name='Raw allele frequency'
  }else{stat_name=stat}
  cov=unique(betai$COVARIABLE_name)
  pop_cov=pop_cov[order(pop_cov[[env_column]]),]
  # Add chr_pos column
  betai$chr_pos=paste(betai$chr, betai$pos, sep="_")
  allele.freq$chr_pos=paste(allele.freq$chr, allele.freq$pos, sep="_")
  # Check there is only results for one covariable in the Betai dataframe
  if(length(cov)>1){stop('Supply Betai dataframe with results from only one covariable')}
  # Select candidates (depending on whether beta is +ve or -ve)
  if(direction=='positive'){
    betai=betai[betai[[beta]]>0,]
  }
  if(direction=='negative'){
    betai=betai[betai[[beta]]<0,]
  }
  if(!(direction %in% c('positive', 'negative', 'both'))){stop("Select beta candidates (either 'positive', 'negative' or 'both')")}
  # If there are any candidates...
  n_cand=nrow(betai)
  if(n_cand>0){
    if(cov%in%c('forest_tree_pct', 'f_over_sum_known_trees')){cov_name='Forest-tree-percentange'}else{cov_name=cov}
    # Select allele frequencies for candidate sites
    allele.freq.select=allele.freq[allele.freq$chr_pos %in% betai$chr_pos,]
    # Add covariable information
    allele.freq.select=merge(allele.freq.select, pop_cov, by='Population')
    # Add Beta values
    allele.freq.select=merge(allele.freq.select, unique(betai[, c('chr_pos', beta)]), all.x=TRUE)
    # Get Beta direction (to colour lines and boxplots)
    allele.freq.select$Direction=NA
    allele.freq.select[allele.freq.select[[beta]]<0, 'Direction']='Negative'
    allele.freq.select[allele.freq.select[[beta]]>0, 'Direction']='Positive'
    # Get correlation coefficient
    allele.freq.select_neg=allele.freq.select[allele.freq.select$Direction=='Negative',]
    allele.freq.select_pos=allele.freq.select[allele.freq.select$Direction=='Positive',]
    if(cor_method=='pearson'){cor_method_name="r"}else{cor_method_name=cor_method}
    subtitle=""
    if(nrow(allele.freq.select_neg)>0){
      cor_res_neg=cor.test(allele.freq.select_neg[[env_column]], 
                           allele.freq.select_neg[[stat]], 
                           method=cor_method)
      subtitle=paste0("r=", signif(cor_res_neg$estimate,2), ", p=",signif(cor_res_neg$p.value,2),"; ")
    }
    if(nrow(allele.freq.select_pos)>0){
      cor_res_pos=cor.test(allele.freq.select_pos[[env_column]], 
                           allele.freq.select_pos[[stat]], 
                           method=cor_method)
      subtitle=paste0(subtitle, "r =", signif(cor_res_pos$estimate,2), ", p=",signif(cor_res_pos$p.value,2))
    }
    if(colour_column!='Direction'){
      allele.freq.select=merge(allele.freq.select, betai[,c("chr_pos", colour_column)], by="chr_pos")
    }
    # Get the average allele frequency for a given SNP in a given habitat (and keep all the relevant columns which are specific to a particular SNP anyway)
    ## This is only really necessary if you happen to have two populations with the same env variable (which of course is far more likely with discrete data)
    allele.freq.select=aggregate(allele.freq.select[[stat]], 
                                 by = list(allele.freq.select$chr_pos, allele.freq.select[[env_column]], allele.freq.select[[beta]], allele.freq.select[[colour_column]]), FUN = mean)
    colnames(allele.freq.select)=c('chr_pos', env_column, beta, colour_column, stat)
    allele.freq.select[[colour_column]]=as.factor(allele.freq.select[[colour_column]])
    
    # Remove NA values
    allele.freq.select=allele.freq.select[!is.na(allele.freq.select[[stat]]),]
    
    # Plot
    cand_plot=ggplot(allele.freq.select, aes_string(x=env_column, y=stat, col=colour_column, fill=colour_column)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=16), 
            plot.subtitle = element_text(hjust = 0.5, size=13),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(title=paste0(title, " (", n_cand, " SNPs)"), 
           subtitle=paste0(subtitle),
           y=paste0(stat_name),
           x=cov_name)
    # Plot line per SNP
    cand_plot=cand_plot+
      geom_line(aes(group=chr_pos), alpha=min(0.5, 150/n_cand)) +
      guides(colour = guide_legend(override.aes = list(alpha=1))) 
    ## If using mean line
    if(big_line=="mean"){
      # Get the average allele frequency of a population
      allele.freq.select_mean=aggregate(allele.freq.select[[stat]], by = list(allele.freq.select$Direction, allele.freq.select[[env_column]]), FUN = mean)
      colnames(allele.freq.select_mean)=c('Direction', 'Population', stat)
      # Plot mean line
      cand_plot=cand_plot+
        geom_line(data=allele.freq.select_mean, aes_string(x='Population', y=stat, group='Direction'), col='black', size=2) +
        geom_line(data=allele.freq.select_mean, aes_string(x='Population', y=stat, col='Direction'), size=1)
    }
    ## If using mean line
    if(big_line=="smooth"){
      cand_plot=cand_plot+
        geom_smooth(method = 'loess', se=F, size=2, col='black') +
        geom_smooth(method = 'loess', se=F, size=1)
    }
    if(big_line=="linear"){
      cand_plot=cand_plot+
        geom_smooth(method='lm', se=F, size=2, col='black') +
        geom_smooth(method='lm', se=F, size=1)
    }
    # If we are adding the position of populations
    if(!is.null(pop_colour_column)){
      pop_cov_top=pop_cov[seq(2, length(pop_cov$Population), 2),]
      pop_cov_bottom=pop_cov[seq(1, length(pop_cov$Population), 2),]
      if(pop_labels){
        cand_plot=cand_plot + 
          scale_x_continuous(breaks = pop_cov_bottom[[env_column]], labels = pop_cov_bottom$Population, 
                             sec.axis = sec_axis(~ ., breaks = pop_cov_top[[env_column]], labels = pop_cov_top$Population))+
          #theme(axis.text.x.top = element_text(angle = 315, vjust = 0, hjust=1, size=7, color = rainbow(nrow(pop_cov_top))),
          #      axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust=1, size=7, color = rainbow(nrow(pop_cov_bottom))))
          theme(axis.text.x.top = element_text(angle = 315, vjust = 0, hjust=1, size=9, color = pop_cov_top[[pop_colour_column]]),
                axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust=1, size=9, color = pop_cov_bottom[[pop_colour_column]],)) + 
          geom_vline(xintercept = pop_cov[[env_column]], col=pop_cov[[pop_colour_column]], linetype='solid', alpha=0.3, size=0.3)
      }else{
        if(stat %in% c("M_P", "DAF", "raw_AF")){min_max=c(-0.04,1.04)}else{min_max=c(min(allele.freq.select$stat), max(allele.freq.select$stat))}
        cand_plot=cand_plot + 
          geom_point(data=pop_cov, aes_string(x=env_column, y='rep(min_max, nrow(pop_cov))[1:nrow(pop_cov)]'), 
                     fill=pop_cov[[pop_colour_column]], 
                     col=pop_cov[[pop_colour_column]],
                     shape=rep(c(24,25), nrow(pop_cov))[1:nrow(pop_cov)], size=3, alpha=1)
      }
    }
    ## If we are dealing with discrete habitats then add boxplots
    if(env_column %in% c('aleman_hab', 'aleman_hab_1cat')){
      cand_plot=cand_plot+geom_boxplot(width=0.3, outlier.size=0.5, outlier.shape=3)
    }
    ## If dealing with raw allele frequencies, set y limits to 0,1 (with some wiggle room so smoothed lines dont dissapear below 0 or above 1)
    if(stat %in% c("M_P", "DAF", "raw_AF")){
      cand_plot=cand_plot + 
        scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
        coord_cartesian(ylim=c(0,1))
    }
    # Legend
    if(colour_column=='Direction'){
      cand_plot=cand_plot+
        scale_discrete_manual(name='Direction', values = c('Negative' = 'red', 'Positive'='deepskyblue'), aesthetics = c("colour")) +
        scale_discrete_manual(name='Direction', values = c('Negative' = 'darkred', 'Positive'='deepskyblue4'), aesthetics = c("fill")) 
    }else{
      cand_plot=cand_plot#+
        #scale_discrete_manual(palette = "Set1", aesthetics = c("colour", "fill"))
        #scale_discrete_manual(values= wes_palette("Darjeeling1", n = length(unique(allele.freq.select[[colour_column]]))), aesthetics = c("colour", "fill")) 
    }
    if(!legend){
      cand_plot=cand_plot + theme(legend.position = "none")
    }
    ## Print plot
    print(cand_plot)
  }else(cat("No sites with ",direction," beta for", unique(betai$COVARIABLE_name), "\n"))
}

run_env_vs_AF=function(betai_cand, stats=c("M_P", "M_Pstd")){
  for(subsp in subsps){
    cat("-", subsp, "\n")
    pop_cov=data.frame(t(cov.in[[subsp]][['Real']]))
    pop_cov$Population=rownames(pop_cov)
    for(stat in stats){
      cat("--", stat, "\n")
      for(fpr in names(betai_cand)){
        for(cov in unique(betai_cand[[fpr]][[subsp]]$COVARIABLE_name)){
          env_vs_AF.v3(betai_cand[[fpr]][[subsp]][betai_cand[[fpr]][[subsp]]$COVARIABLE_name==cov,], 
                       allele.freq[[subsp]], 
                       pop_cov, 
                       stat=stat, 
                       direction='both', 
                       subtitle=paste0("Candidates: FPR=", fpr), 
                       beta='M_Beta.median')
        }
      }
    }
  }
}

run_env_vs_AF.v2=function(betai, stats=c("M_P", "M_Pstd"), fprs=c(0.005, 0.0005, 0.00005)){
  for(subsp in subsps){
    cat("-", subsp, "\n")
    pop_cov=data.frame(t(cov.in[[subsp]][['Real']]))
    pop_cov$Population=rownames(pop_cov)
    for(stat in stats){
      cat("--", stat, "\n")
      for(fpr in fprs){
        for(cov in unique(betai[[subsp]]$COVARIABLE_name)){
          env_vs_AF.v3(betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr,], 
                       allele.freq[[subsp]], 
                       pop_cov, 
                       stat=stat, 
                       direction='both', 
                       subtitle=paste0("Candidates: FPR=", fpr), 
                       beta='M_Beta.median')
        }
      }
    }
  }
}

run_env_vs_AF.v4=function(betai_list, cov.in, allele.freq, stats=c("M_P", "M_Pstd"), fprs=c(0.005, 0.001, 0.0005), 
                          big_line="smooth", colour_column='Direction', legend=TRUE, pop_labels=TRUE){
  subsps=names(betai_list)
  for(stat in stats){
    for(subsp in subsps){
      #cat("-", subsp, "\n")
      if(subsp=="all"){subsp_name="All subspecies"}
      if(subsp=="c"){subsp_name="Central"}
      if(subsp=="e"){subsp_name="Eastern"}
      if(subsp=="ce"){subsp_name="Central-Eastern"}
      if(subsp=="n"){subsp_name="Nigeria-Cameroon"}
      if(subsp=="w"){subsp_name="Western"}
      pop_cov=data.frame(t(cov.in[[subsp]][['Real']]))
      pop_cov$Population=rownames(pop_cov)
      pop_cov=pop_cov[order(pop_cov$Population),]
      pops=unique(pop_cov$Population)
      pops=pops[order(pops)]
      pop_cov$subsp_col=NA
      for(pop in pop_cov$Population){
        if(startsWith(pop, "c.")){pop_cov[pop_cov$Population == pop, 'subsp_col']="green2"}
        if(startsWith(pop, "e.")){pop_cov[pop_cov$Population == pop, 'subsp_col']="darkorange"}
        if(startsWith(pop, "n.")){pop_cov[pop_cov$Population == pop, 'subsp_col']="red3"}
        if(startsWith(pop, "w.")){pop_cov[pop_cov$Population == pop, 'subsp_col']="blue"}
      }
      #cat("--", stat, "\n")
      for(fpr in fprs){
        for(cov in unique(betai_list[[subsp]]$COVARIABLE_name)){
          env_vs_AF.v4(betai_list[[subsp]][betai_list[[subsp]]$COVARIABLE_name==cov & betai_list[[subsp]]$fpr<=fpr,], 
                       allele.freq[[subsp]], 
                       pop_cov, 
                       env_column=cov,
                       pop_colour_column='subsp_col',
                       stat=stat, 
                       direction='both', 
                       title=paste0(subsp_name,", FPR<", 100*fpr, "%"), 
                       beta='M_Beta.median',
                       big_line=big_line,
                       colour_column=colour_column,
                       legend=legend,
                       pop_labels=pop_labels)
        }
      }
    }
  }
}

env_AF_diff=function(betai, allele.freq, pop_cov, stat="M_Pstd", direction='both', subtitle='', beta='M_Beta', largest_diff=FALSE){
  cov=unique(betai$COVARIABLE_name)
  # Checks
  ## Check there is only results for one covariable in the Betai dataframe
  if(length(cov)>1){stop('Supply Betai dataframe with results from only one covariable')}
  ## Make a file with just a row per SNP (no gene annotation)
  betai=unique(betai[colnames(betai)!='gene'])
  if(nrow(betai) != length(unique(betai$MRK))){stop("Betai dataframe contains multiple rows corresponding to single SNPs")}
  # Select candidates (depending on whether beta is +ve or -ve)
  if(direction=='positive'){betai=betai[betai[[beta]]>0,]}
  if(direction=='negative'){betai=betai[betai[[beta]]<0,]}
  if(!(direction %in% c('positive', 'negative', 'both'))){stop("Select beta candidates (either 'positive', 'negative' or 'both')")}
  # If there are any SNPs...
  if(nrow(betai)>0){
    # Select allele frequencies for candidate sites
    allele.freq.select=allele.freq[allele.freq$MRK %in% betai$MRK,]
    # Add covariable information
    allele.freq.select=merge(allele.freq.select, pop_cov, by='Population')
    # Add Beta values
    allele.freq.select=merge(allele.freq.select, unique(betai[, c('MRK', beta)]), all.x=TRUE)
    # Get Beta direction (to colour lines and boxplots)
    allele.freq.select$Direction=NA
    allele.freq.select[allele.freq.select[[beta]]<0, 'Direction']='Negative'
    allele.freq.select[allele.freq.select[[beta]]>0, 'Direction']='Positive'
    # Change non-focal habitat labels to 'other'
    allele.freq.select[allele.freq.select$Vrishti_habitats!=cov, 'Vrishti_habitats']='Other'
    # Focus only on the largest difference between focal habitat and others 
    if(largest_diff==TRUE){
      ## Change the allele frequency column name to AF (just for this section to make it easier to use dplyr)
      colnames(allele.freq.select)[colnames(allele.freq.select)==stat]="AF"
      ## Select SNPs with a NEGATIVE association
      allele.freq.select.neg=allele.freq.select[allele.freq.select$Direction=='Negative',]
      ### Get the MINIMUM allele frequency of any population from the FOCAL habitat
      allele.freq.select.neg.min=allele.freq.select.neg[allele.freq.select.neg$Vrishti_habitats==cov,] %>% 
        group_by(MRK, Direction, Vrishti_habitats) %>% 
        summarise(AF = min(AF))
      ### Get the MAXIMUM allele frequency of any population from the OTHER habitat
      allele.freq.select.neg.max=allele.freq.select.neg[allele.freq.select.neg$Vrishti_habitats=='Other',] %>% 
        group_by(MRK, Direction, Vrishti_habitats) %>% 
        summarise(AF = max(AF))
      ## Select SNPs with a POSITIVE association
      allele.freq.select.pos=allele.freq.select[allele.freq.select$Direction=='Positive',]
      ### Get the MINIMUM allele frequency of any population from the OTHER habitat
      allele.freq.select.pos.min=allele.freq.select.pos[allele.freq.select.pos$Vrishti_habitats=='Other',] %>% 
        group_by(MRK, Direction, Vrishti_habitats) %>% 
        summarise(AF = min(AF))
      ### Get the MAXIMUM allele frequency of any population from the FOCAL habitat
      allele.freq.select.pos.max=allele.freq.select.pos[allele.freq.select.pos$Vrishti_habitats==cov,] %>% 
        group_by(MRK, Direction, Vrishti_habitats) %>% 
        summarise(AF = max(AF))
      ## Add all these datasets back together
      ### We now have a row per SNP with allele frequencies from the two populations with the biggest difference for that SNP
      allele.freq.select=rbind(allele.freq.select.neg.min, allele.freq.select.neg.max, allele.freq.select.pos.min, allele.freq.select.pos.max)
      ### Change the allele frequency column name back
      colnames(allele.freq.select)[colnames(allele.freq.select)=="AF"]=paste0(stat)
    }
    # Get the average allele frequency for a given SNP in a given habitat (and keep all the relevant columns which are specific to a particular SNP anyway)
    ## This does nothing if largest_diff==TRUE
    allele.freq.select=aggregate(allele.freq.select[[stat]], 
                                 by = list(allele.freq.select$MRK, allele.freq.select$Vrishti_habitats, allele.freq.select$Direction), FUN = mean)
    colnames(allele.freq.select)=c('MRK', 'Vrishti_habitats', 'Direction', stat)
    # Get difference in allele frequencies between focal habitat type and all others 
    allele.freq.select=reshape(allele.freq.select, idvar = c("MRK", "Direction"), timevar = "Vrishti_habitats", direction = "wide")
    allele.freq.select$Difference=allele.freq.select[[paste0(stat,".",cov)]]-allele.freq.select[[paste0(stat,".Other")]]
    # Plot
    cand_plot=ggplot(allele.freq.select, aes(x=Difference, col=Direction, fill=Direction)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=25), plot.subtitle = element_text(hjust = 0.5, size=15)) +
      labs(title=paste0(cov, " Candidates: ", stat), 
           subtitle=paste0(subtitle, " (", nrow(betai), " SNPs)"),
           x=paste0("Mean ", stat, " in ", cov, " Populations - Mean ", stat, " in Other Populations"),
           y="Density") +
      geom_density(alpha=0.3) +
      scale_discrete_manual(name='Direction', values = c('Negative' = 'red', 'Positive'='green3'), aesthetics = c("colour", "fill")) 
    print(cand_plot)
  }else(cat("No sites with ",direction," beta for", unique(betai$COVARIABLE_name), "\n"))
}

aleman_hab_maps=function(env_file){
  # Load maps
  world_map=map_data("world")
  chimp.ranges=readOGR("../../../maps/chimp_ranges/chimp_ranges.shp")
  aleman_hab_map=raster("../../../environmental_data/input/aleman_et.el.2020/biome_potential.tif")
  #aleman_hab_raster.legend=data.frame(biome=c('savanna', 'savanna_mosaic', 'forest_mosaic', 'forest'),
  #                  code=1:4)
  aleman_hab_map=as(aleman_hab_map, "SpatialPixelsDataFrame")
  aleman_hab_map=as.data.frame(aleman_hab_map)
  colnames(aleman_hab_map)=c("Aleman et al. Habitats", "lat", "long")
  aleman_hab_map[aleman_hab_map$`Aleman et al. Habitats`==4, 'Aleman et al. Habitats']='forest'
  aleman_hab_map[aleman_hab_map$`Aleman et al. Habitats`==3, 'Aleman et al. Habitats']='mosaic_forest'
  aleman_hab_map[aleman_hab_map$`Aleman et al. Habitats`==2, 'Aleman et al. Habitats']='mosaic_savanna'
  aleman_hab_map[aleman_hab_map$`Aleman et al. Habitats`==1, 'Aleman et al. Habitats']='savanna'
  # Limits
  ## x
  xmin=min(env_file$Longitude)
  xmax=max(env_file$Longitude)
  xmin=xmin-(xmax-xmin)*0.1
  xmax=xmax+(xmax-xmin)*0.1
  ## x
  ymin=min(env_file$Latitude)
  ymax=max(env_file$Latitude)
  ymin=ymin-(ymax-ymin)*0.15
  ymax=ymax+(ymax-ymin)*0.15
  
  # Raster plot
  raster=ggplot(NULL) +
    # Map
    theme_void() +
    #geom_polygon(data=world_map, aes(x = long, y = lat, group = group), fill="grey93", colour = "darkgrey") +
    geom_tile(data=aleman_hab_map, aes(x = lat, y = long, fill = `Aleman et al. Habitats`), alpha=0.9) +
    #coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax))+
    scale_discrete_manual(name="Aleman et al.,\n(2020)", 
                          values = c('forest' = 'darkgreen', 'mosaic_forest'='yellow4', 'mosaic_savanna' = 'orange3', 'savanna'='goldenrod1'), 
                          aesthetics = c("colour", "fill"))
  print(raster)
  raster_ranges=raster +
    # Subspecies ranges
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="troglodytes",], aes(x = long, y = lat, group = group), col="white", alpha=0, size=1.5) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="troglodytes",], aes(x = long, y = lat, group = group), col="green3",alpha=0) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="schweinfurthii",], aes(x = long, y = lat, group = group), col="white",alpha=0, size=1.5) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="schweinfurthii",], aes(x = long, y = lat, group = group), col="darkorange",alpha=0) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="ellioti",], aes(x = long, y = lat, group = group), col="white",alpha=0, size=1.5) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="ellioti",], aes(x = long, y = lat, group = group), col="red",alpha=0) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="verus",], aes(x = long, y = lat, group = group), col="white", alpha=0, size=1.5) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="verus",], aes(x = long, y = lat, group = group), col="blue", alpha=0)
  print(raster_ranges)
  raster_ranges_points=raster_ranges +
    # Plot points
    geom_point(data=env_file, aes(x=Longitude, y=Latitude), col='white', size=3, shape=15) +
    geom_point(data=env_file, aes(x=Longitude, y=Latitude), col='black', size=2, shape=15) 
  print(raster_ranges_points)
  # Point plot
  ggplot(NULL) +
    # Map
    theme_void() +
    geom_polygon(data=world_map, aes(x = long, y = lat, group = group), fill="grey93", colour = "darkgrey") +
    #coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax))+
    coord_sf(xlim=layer_scales(raster_ranges_points)$x$range$range, ylim=layer_scales(raster_ranges_points)$y$range$range)+
    # Subspecies ranges
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="troglodytes",], aes(x = long, y = lat, group = group), fill="green3", alpha=0.3) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="schweinfurthii",], aes(x = long, y = lat, group = group), fill="darkorange", alpha=0.3) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="ellioti",], aes(x = long, y = lat, group = group), fill="red", alpha=0.3) +
    geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="verus",], aes(x = long, y = lat, group = group), fill="blue", alpha=0.3) +
    # Plot points
    geom_point(data=env_file, aes(x=Longitude, y=Latitude, col=aleman_hab), size=3, shape=15) +
    scale_discrete_manual(name="Aleman et al.,\n(2020)", 
                          values = c('forest' = 'darkgreen', 'mosaic_forest'='yellow4', 'mosaic_savanna' = 'orange3', 'savanna'='goldenrod1'), 
                          aesthetics = c("colour", "fill"))
}

aleman_hab_maps_random=function(env_file, cov.in){
  # Load maps
  world_map=map_data("world")
  chimp.ranges=readOGR("../../../maps/chimp_ranges/chimp_ranges.shp")
  # Limits
  ## x
  xmin=min(env_file$Longitude)
  xmax=max(env_file$Longitude)
  xmin=xmin-(xmax-xmin)*0.1
  xmax=xmax+(xmax-xmin)*0.1
  ## x
  ymin=min(env_file$Latitude)
  ymax=max(env_file$Latitude)
  ymin=ymin-(ymax-ymin)*0.15
  ymax=ymax+(ymax-ymin)*0.15
  
  # No need to loop over all subspecies datasets as they use the same random habitat assignments as plotted here
  ## i.e. for the random 1 repeat for western, I just selected the habitats for the western subspecies
  for(subsp in 'all'){
    for(rep in names(cov.in[[subsp]])){
      # Reformat covariable file
      hab.df=data.frame(t(cov.in[[subsp]][[rep]]))
      hab.df$Habitat=NA
      hab.df[hab.df$forest==1, 'Habitat']='forest'
      hab.df[hab.df$mosaic_forest==1, 'Habitat']='mosaic_forest'
      hab.df[hab.df$mosaic_savanna==1, 'Habitat']='mosaic_savanna'
      hab.df[hab.df$savanna==1, 'Habitat']='savanna'
      ## Add coordinates
      hab.df=merge(unique(env_file[,c('Population', 'Longitude', 'Latitude')]), hab.df, by.x='Population', by.y=0)
      # Plot
      map=ggplot(NULL) +
        ## Map
        theme_void() +
        geom_polygon(data=world_map, aes(x = long, y = lat, group = group), fill="grey93", colour = "darkgrey") +
        coord_sf(xlim=c(xmin, xmax), ylim=c(ymin, ymax)) +
        labs(title=paste0(rep)) +
        theme(plot.title = element_text(hjust = 0.5, size=25)) +
        ## Subspecies ranges
        geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="troglodytes",], aes(x = long, y = lat, group = group), fill="green3", alpha=0.3) +
        geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="schweinfurthii",], aes(x = long, y = lat, group = group), fill="darkorange", alpha=0.3) +
        geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="ellioti",], aes(x = long, y = lat, group = group), fill="red", alpha=0.3) +
        geom_polygon(data=chimp.ranges[chimp.ranges$subspecies=="verus",], aes(x = long, y = lat, group = group), fill="blue", alpha=0.3) +
        ## Plot points
        geom_point(data=hab.df, aes(x=Longitude, y=Latitude, col=Habitat), size=3, shape=15) +
        scale_discrete_manual(name="Habitats", 
                              values = habitat_col, aesthetics = c("colour", "fill"))
      print(map)
    }
  }
}

env_var_map=function(df, title=NULL, col=NULL){
  if(!all(c("Population", "Latitude", "Longitude") %in% colnames(df))){
    stop("Input must have at least four columns ; 'Population', 'Latitude', 'Longitude' and then a column for each environmnetal variable to be plotted")
  }
  # Load maps
  world_map=map_data("world")
  chimp.ranges=readOGR("~/OneDrive - University College London/Projects/PanAf/maps/chimp_ranges/chimp_ranges.shp")
  ## Downloaded version 5.0.0 from https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-rivers-lake-centerlines/
  rivers=readOGR("~/OneDrive - University College London/Projects/PanAf/maps/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp")
  rivers=rivers[grepl("Congo|Niger|Ubangi|Ogooué|Tanganyika|Sanaga|Benue|Uelé|Lualaba|Chambeshi|Chinko|Mbomou|Uele", rivers$name_en),]
  ## Download water bodies Zip (last updated Nov 14, 2018) https://datacatalog.worldbank.org/search/dataset/0040797
  lakes=readOGR("~/OneDrive - University College London/Projects/PanAf/maps/africawaterbody/Africa_waterbody.shp")
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
  
  for(env_var in colnames(df[, !colnames(df) %in% c("Population", "Latitude", "Longitude"), drop=FALSE])){
    # name if plotting forest tree percentage 
    if(env_var=="forest_tree_pct"){legend_name="Forest\ntrees (%)"}else{legend_name=env_var}
    # Plot
    map=ggplot(NULL) +
      ## Map
      theme_void() +
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
      ## Population labels
      geom_text_repel(data=df, aes(x=Longitude, y=Latitude, label=Population), 
                      point.padding=1,
                      #fontface= "bold",
                      min.segment.length = unit(0, 'lines'),
                      size=2.5,
                      seed=1,
                      max.overlaps = Inf,
                      bg.color = "white",
                      bg.r = 0.1,
                      fontface = 'bold') +
      ## Plot points
      geom_point(data=df, aes_string(x='Longitude', y='Latitude'), col='black', size=6, shape=16) +
      geom_point(data=df, aes_string(x='Longitude', y='Latitude', col=env_var), size=5, shape=16) +
      ## Population labels (once more but with no lines so the text overlays the points but lines do not)
      geom_text_repel(data=df, aes(x=Longitude, y=Latitude, label=Population), 
                      point.padding=1,
                      #fontface= "bold",
                      min.segment.length = unit(999, 'lines'),
                      size=2.5,
                      seed=1,
                      max.overlaps = Inf,
                      bg.color = "white",
                      bg.r = 0.1,
                      fontface = 'bold')
    
    if(is.numeric(df[[env_var]])){
      if(is.null(col)){colour=viridis(100)}else{colour=col}
      map=map+scale_colour_gradientn(name=legend_name, colours = colour, aesthetics = c("colour", "fill"))
    }else{
      if(is.null(col)){colour=viridis(length(unique(df[[env_var]])))}else{colour=col}
      map=map+scale_discrete_manual(name=legend_name, values = colour, aesthetics = c("colour", "fill"))
    }
    print(map)
  }
}

plot_Pdelta=function(Pdelta){
  subsps=names(Pdelta)
  for(subsp in subsps){
    for(stat in c('M_P.median')){
      ## Select only results from real (not random) runs
      Pdelta_plot_df=Pdelta[[subsp]][Pdelta[[subsp]]$Run=='real',]
      ## Get mean across subsets
      Pdelta_plot_df=aggregate(Pdelta_plot_df[stat], Pdelta_plot_df[c("COVARIABLE_name")], mean)
      plot=ggplot(Pdelta_plot_df, aes_string(x='COVARIABLE_name', y=stat, col='COVARIABLE_name', fill='COVARIABLE_name')) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=25),
              axis.text.x = element_text(angle = 45, hjust=1)) +
        labs(title=paste0(subsp, " ", stat)) +
        geom_bar(stat="identity", col='black') +
        scale_discrete_manual(name='Habitat', values = habitat_col, aesthetics = c("colour", "fill"))
      # Extract legend
      legend=get_legend(plot)
      # If there are too many covariables, plot the legend separately
      if(length(covs)>10){
        plot=plot+theme(legend.position = "none")
      }
      print(plot)
    }
  }
  if(length(covs)>10){
    grid.newpage()
    grid.draw(legend)
  }
}

plot_Pdelta.subsp_panels=function(Pdelta){
  subsps=names(Pdelta)
  Pdelta_plot_df=NULL
  for(stat in c('M_P.median')){
    for(subsp in subsps){
      Pdelta[[subsp]]$Subspecies=subsp
      ## Select only results from real (not random) runs
      Pdelta_plot_df.subsp=Pdelta[[subsp]][Pdelta[[subsp]]$Run=='real',]
      if(!is.null(Pdelta_plot_df)){
        Pdelta_plot_df=rbind(Pdelta_plot_df, Pdelta_plot_df.subsp)
      }else{
        Pdelta_plot_df=Pdelta_plot_df.subsp
      }
    }
    ## Get mean across subsets
    Pdelta_plot_df=aggregate(Pdelta_plot_df[stat], Pdelta_plot_df[c("Subspecies", "COVARIABLE_name")], mean) 
      plot=ggplot(Pdelta_plot_df, aes_string(x='COVARIABLE_name', y=stat, col='COVARIABLE_name', fill='COVARIABLE_name')) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=25),
              axis.text.x = element_text(angle = 45, hjust=1)) +
        labs(title=paste0(subsp, " ", stat)) +
        geom_bar(stat="identity", col='black') +
        facet_wrap(~Subspecies) +
        scale_discrete_manual(name='Habitat', values = habitat_col, aesthetics = c("colour", "fill"))
      # Extract legend
      legend=get_legend(plot)
      # If there are too many covariables, plot the legend separately
      if(length(covs)>10){
        plot=plot+theme(legend.position = "none")
      }
        print(plot)
  }
  if(length(covs)>10){
    grid.newpage()
    grid.draw(legend)
  }
}

plot_Pdelta_random=function(Pdelta){
  subsps=names(Pdelta)
  for(subsp in subsps){
    Pdelta_plot_df=Pdelta[[subsp]]
    # Add run type column
    Pdelta_plot_df$run_type='Real'
    Pdelta_plot_df[Pdelta_plot_df$Run!='real', 'run_type']='Random'
    # Get mean across subsets
    Pdelta_plot_df=aggregate(Pdelta_plot_df['M_P.median'], Pdelta_plot_df[c("Run", "COVARIABLE_name", "run_type")], mean)
    #print(Pdelta[[subsp]][order(Pdelta[[subsp]]$COVARIABLE_name), c('COVARIABLE_name', 'Run', 'M_P.median')])
    runs=unique(Pdelta_plot_df$Run)
    runs=runs[order(runs)]
    Pdelta_plot_df$Run=factor(Pdelta_plot_df$Run, levels=c("real", runs[runs!="real"]))
    for(hab in unique(Pdelta_plot_df$COVARIABLE_name)){
      M_P_plot=ggplot(Pdelta_plot_df[Pdelta_plot_df$COVARIABLE_name==hab,], aes(x=Run, y=M_P.median, fill=COVARIABLE_name, col=COVARIABLE_name, pattern=run_type)) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="none") +
        geom_bar_pattern(stat="identity",
                         col = "black", 
                         pattern_fill = "black",
                         pattern_angle = 45,
                         pattern_density = 0.1,
                         pattern_spacing = 0.025,
                         pattern_key_scale_factor = 0.6) +
        geom_hline(yintercept=mean(Pdelta_plot_df[Pdelta_plot_df$COVARIABLE_name==hab & Pdelta_plot_df$Run!='real', 'M_P.median']), col="black", linetype="twodash", size=1, alpha=0.75) +
        labs(title=paste0(subsp, " ", hab, " : Median M_P"),
             subtitle = "Dashed Line Shows Mean Across Random Median Values",) +
        scale_discrete_manual("", values = habitat_col, aesthetics = c("colour", "fill")) +
        scale_pattern_manual(values = c(Real = "none", Random = "stripe"))
      print(M_P_plot)
    }
  }
}

plot_Pdelta_chr21=function(Pdelta, Pdelta_chr21){
  subsps=names(Pdelta)
  for(subsp in subsps){
    Pdelta_plot_df_ex=Pdelta[[subsp]]
    Pdelta_plot_df_chr21=Pdelta_chr21[[subsp]]
    # Add run type column
    ## Exome
    ### Select only real data
    Pdelta_plot_df_ex=Pdelta_plot_df_ex[Pdelta_plot_df_ex$Run=='real',]
    Pdelta_plot_df_ex$Dataset='Exome'
    ### Remove run column
    Pdelta_plot_df_ex=Pdelta_plot_df_ex[,'Run' != colnames(Pdelta_plot_df_ex)]
    ## chr21
    Pdelta_plot_df_chr21$Dataset='chr21'
    # Combine datasets
    Pdelta_plot_df=rbind(Pdelta_plot_df_ex, Pdelta_plot_df_chr21)
    # Get mean across subsets
    Pdelta_plot_df=aggregate(Pdelta_plot_df['M_P.median'], Pdelta_plot_df[c("COVARIABLE_name", "Dataset")], mean)
    Pdelta_plot_df$Dataset=factor(Pdelta_plot_df$Dataset, levels=c("Exome", "chr21"))
    for(hab in unique(Pdelta_plot_df$COVARIABLE_name)){
      M_P_plot=ggplot(Pdelta_plot_df[Pdelta_plot_df$COVARIABLE_name==hab,], aes(x=Dataset, y=M_P.median, fill=COVARIABLE_name, col=COVARIABLE_name, pattern=Dataset)) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="none") +
        geom_bar_pattern(stat="identity",
                         col = "black", 
                         pattern_fill = "black",
                         pattern_angle = 45,
                         pattern_density = 0.1,
                         pattern_spacing = 0.025,
                         pattern_key_scale_factor = 0.6) +
        labs(title=paste0(subsp, " ", hab, " : Median M_P")) +
        scale_discrete_manual("", values = habitat_col, aesthetics = c("colour", "fill")) +
        scale_pattern_manual(values = c(Exome = "none", chr21 = "stripe"))
      print(M_P_plot)
    }
  }
}

plot_betai_stats=function(betai, cols=habitat_col){
  subsp=names(betai)
  # Plot distributions of all stats
  for(subsp in subsps){
    cat(subsp)
    for(stat in c('M_Beta', 'SD_Beta', 'PIP', '`BF(dB)`')){
      plot=ggplot(betai[[subsp]], aes_string(x=stat, col='COVARIABLE_name', fill='COVARIABLE_name')) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=25)) +
        geom_density(alpha=0.1) +
        scale_discrete_manual(name='', values = cols, aesthetics = c("colour", "fill")) +
        labs(title=paste0(subsp, " ", stat), y='Density')
      ## If plotting beta, use a (pseudo) log y axis
      if(stat=='M_Beta' | stat=='PIP'){
        plot=plot+scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), expand = c(0, 0)) + labs(y=expression(log[10]~density))
      }
      # Extract legend
      legend=get_legend(plot)
      # If there are too many covariables, plot the legend separately
      if(length(covs)>10){
        plot=plot+theme(legend.position = "none")
      }
      print(plot)
    }
  }
  if(length(covs)>10){
    grid.newpage()
    grid.draw(legend)
  }
}

plot_betai_stats.v2=function(betai, null=NULL, cols=habitat_col, stats=c('M_Beta.median', '`BF(dB).median`')){
  subsp=names(betai)
  # Plot distributions of all stats
  for(subsp in subsps){
    cat(subsp)
    for(stat in stats){
      plot=ggplot(NULL) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=25)) +
        geom_density(data=betai[[subsp]], aes_string(x=stat, col='COVARIABLE_name', fill='COVARIABLE_name', linetype = shQuote("Exome")), alpha=0.1) +
        labs(title=paste0(subsp, " ", stat), y='Density') 
      ## If plotting beta, use a (pseudo) log y axis
      if(grepl('M_Beta', stat, fixed = TRUE) | grepl('PIP', stat, fixed = TRUE)){
        plot=plot+scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), expand = c(0, 0)) + labs(y=expression(log[10]~(Density)))
      }
      # Extract legend
      legend=get_legend(plot)
      # If there are too many covariables, plot the legend separately
      if(length(covs)>10){
        plot=plot+theme(legend.position = "none")
      }
      # Add null distribution?
      if(!is.null(null)){
        plot=plot+
          geom_density(data=null[[subsp]], aes_string(x=stat, col='COVARIABLE_name', fill='COVARIABLE_name', linetype = shQuote("Null")), alpha=0.1) +
          scale_linetype_manual(name='', values = c("Exome"="solid", "Null"="dashed"))
      }
      plot=plot+
        scale_discrete_manual(name='', values = c(cols), aesthetics = c("colour", "fill"))
      print(plot)
    }
  }
  if(length(covs)>10){
    grid.newpage()
    grid.draw(legend)
  }
}

core_vs_aux=function(betai){
  subsps=names(betai)
  for(subsp in subsps){
    cat(subsp, "\n")
    for(cov in unique(betai[[subsp]]$COVARIABLE_name)){
      for(stat in c('M_Beta.median', '`abs(M_Beta.median)`', '`BF(dB).median`')){
        plot_data=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov,]
        if(stat=='`abs(M_Beta.median)`'){
          plot_data$`abs(M_Beta.median)`=abs(plot_data$M_Beta.median)
        }
        plot=ggplot(data=plot_data, aes_string(x="XtXst.med", y=stat)) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size=25)) +
          geom_hex(bins = 50) +
          geom_smooth(col='black') +
          scale_fill_gradientn(colours = rev(rainbow(5)), trans = "log") +
          labs(title=paste0(subsp, " ", cov, ":\n", stat))
        print(plot)
      }
    }
  }
}

bf_thresh_cand=function(betai, bf.thresh, habitat_col){
  subsps=names(betai)
  for(subsp in subsps){
    cand_freq=data.frame(table(betai[[subsp]][as.numeric(betai[[subsp]]$`BF(dB).median`)>bf.thresh, 'COVARIABLE_name']))
    colnames(cand_freq)=c("Habitat", "Number of Candidates")
    cand_freq_plot=ggplot(cand_freq, aes(x=Habitat, y=`Number of Candidates`, fill=Habitat))+
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust=1)) +
      geom_bar(stat="identity", col='black')+
      labs(title=paste0(subsp, ": BF > ", bf.thresh)) +
      scale_discrete_manual(name="Habitats", 
                            values = habitat_col, 
                            aesthetics = c("colour", "fill"))
    # If there are lost of categories then remove the legend
    if(nrow(cand_freq>10)){
      cand_freq_plot=cand_freq_plot + theme(legend.position = "none")
    }
    print(cand_freq_plot)
  }
}

bf_thresh_cand_subsp_panels=function(betai, bf.thresh, habitat_col){
  subsps=names(betai)
  # Combine subspecies datasets
  cand_freq=NULL
  for(subsp in subsps){
    betai[[subsp]]$chr_pos=paste(betai[[subsp]]$chr, betai[[subsp]]$pos, sep="_")
    betai[[subsp]]=unique(betai[[subsp]][c("chr_pos", "BF(dB).median", "COVARIABLE_name")])
    cand_freq.subsp=data.frame(table(betai[[subsp]][as.numeric(betai[[subsp]]$`BF(dB).median`)>bf.thresh, 'COVARIABLE_name']))
    cand_freq.subsp$Subspecies=subsp
    cand_freq.subsp$Tail_pct=100*nrow(betai[[subsp]][as.numeric(betai[[subsp]]$`BF(dB).median`)>bf.thresh, ])/nrow(betai[[subsp]])
    if(!is.null(cand_freq)){
      cand_freq=rbind(cand_freq, cand_freq.subsp)
    }else{
      cand_freq=cand_freq.subsp
    }
  }
  print(cand_freq)
  # Plot
  colnames(cand_freq)=c("Habitat", "Number of Candidates", "Subspecies")
  cand_freq_plot=ggplot(cand_freq, aes(x=Habitat, y=`Number of Candidates`, fill=Habitat))+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust=1)) +
    geom_bar(stat="identity", col='black')+
    labs(title=paste0("BF > ", bf.thresh)) +
    scale_discrete_manual(name="Habitats", 
                          values = habitat_col, 
                          aesthetics = c("colour", "fill")) +
    facet_wrap(~Subspecies)
  # If there are lost of categories then remove the legend
  if(nrow(cand_freq>10)){
    cand_freq_plot=cand_freq_plot + theme(legend.position = "none")
  }
  print(cand_freq_plot)
}

combine_random_runs=function(betai, random_reps, vol_plots=FALSE){
  subsps=names(betai)
  betai.rand.total=list()
  for(subsp in subsps){
    cat("-", subsp, "\n-- Random Runs\n")
    betai.rand.total[[subsp]]=list()
    for(rep in random_reps[[subsp]]){
      # Plot volcano plot
      if(vol_plots){
        aux_volcano.v2(betai[[subsp]], x=paste0('`M_Beta.median_RANDOM_rep', rep,"`"), y=paste0('`BF(dB).median_RANDOM_rep', rep,"`"), colour_col='COVARIABLE_name', colours=habitat_col)
      }
      ## Extract relevant columns
      betai.rand.total[[subsp]][[rep]]=betai[[subsp]][grepl(paste0("^chr$|^pos$|^COVARIABLE_name$|^coverage$|.*median_RANDOM_rep", rep), colnames(betai[[subsp]]))]
      colnames(betai.rand.total[[subsp]][[rep]])=gsub(paste0("_RANDOM_rep",rep), "", colnames(betai.rand.total[[subsp]][[rep]]))
    }
    cat("-- Combined Random Runs\n")
    # Stack data frame from random runs 
    betai.rand.total[[subsp]]=do.call(rbind, betai.rand.total[[subsp]])
    # Plot all repeats combined (this is used to define FPR later)
    if(vol_plots){
      aux_volcano.v2(betai.rand.total[[subsp]], x='M_Beta.median', y='`BF(dB).median`', colour_col='COVARIABLE_name', colours=habitat_col)
    }
    plot=ggplot(betai.rand.total[[subsp]], aes(x=`BF(dB).median`, col=COVARIABLE_name, fill=COVARIABLE_name)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size=25)) +
      geom_density(alpha=0.1) +
      scale_discrete_manual(name='Habitat', values = habitat_col, aesthetics = c("colour", "fill")) +
      labs(title=paste0(subsp, ": null BF distribution"), y='Density')
    print(plot)
  }
  return(betai.rand.total)
}

select_candidates=function(betai, betai.rand.total, fprs=c(0.005, 0.0005, 0.00005), n_bins=15, top.cov=0.025, habitat_col){
  betai_cand=list()
  cover_bins=list()
  subsps=names(betai)
  for(subsp in subsps){
    # Assign SNPs to coverage bins
    
    # Get all coverage values
    cover=betai[[subsp]]$coverage
    # Ensure that the lowest and highest bins contain 5% of the data points 
    low=max(head(cover[order(cover)], 0.05*length(cover)))
    high=min(head(cover[order(-cover)], 0.05*length(cover)))
    ## Separate into n bins between the low and high bin thresholds so we end up with n+2 bins over all (when you include the highest and lowest bins)
    cover_bins[[subsp]]=c(0, seq(from = low, to = high, by = (high-low)/(n_bins-2)), max(cover))
    # Assign each SNP to a bin
    betai[[subsp]]$coverage_bin=NA
    for(i in 1:(length(cover_bins[[subsp]])-1)){
      betai[[subsp]][betai[[subsp]]$coverage > cover_bins[[subsp]][i] & betai[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
    }
  # Get BF Threshold per Bin
  for(fpr in fprs){
    fpr_pct=paste0(100*fpr, "pct")
    # Make chr_pos column
    betai[[subsp]]$chr_pos=paste(betai[[subsp]]$chr, betai[[subsp]]$pos, sep="_")
    betai.rand.total[[subsp]]$chr_pos=paste(betai.rand.total[[subsp]]$chr, betai.rand.total[[subsp]]$pos, sep="_")
    # Get list of habitat names
    habs=unique(betai[[subsp]]$COVARIABLE_name)
    habs=habs[order(habs)]
    # Prime the column for BF threshold per bin per habitat type for a particular FPR
    betai[[subsp]]$tmp=0
    colnames(betai[[subsp]])[colnames(betai[[subsp]])=='tmp']=paste0(fpr_pct, '_FPR_BF_thresh')
    # Loop over each bin
    bins=unique(betai[[subsp]]$coverage_bin)
    for(bin in bins[order(bins)]){
      ## Get the chr_pos for SNPs in the bin
      bin_snp=betai[[subsp]][betai[[subsp]]$coverage_bin == bin,]$chr_pos
      null_BFs=c()
      ### Extract SNPs from the bin in the random runs (and combine the information across random runs)
      for(hab in habs){
        #### Get null BFs by combing all BFs from SNPs in the coverage bin for the given habitat type
        null_BFs=betai.rand.total[[subsp]][(betai.rand.total[[subsp]]$chr_pos %in% bin_snp)
                                           & (betai.rand.total[[subsp]]$COVARIABLE_name==hab), 'BF(dB).median']
        ## Get BF threshold for given fpr
        bf.thresh=min(head(null_BFs[order(-null_BFs)], length(null_BFs)*fpr))
        betai[[subsp]][betai[[subsp]]$coverage_bin==bin & betai[[subsp]]$COVARIABLE_name==hab, paste0(fpr_pct, '_FPR_BF_thresh')]=bf.thresh
        ## Report thresholds per bin
        #cat("bin ", bin, " ", hab, " BF thresh: ", bf.thresh, "\n")
      }
    }
    # Plot thresholds
    thresh_df=unique(betai[[subsp]][c('coverage_bin', 'COVARIABLE_name', paste0(fpr_pct, '_FPR_BF_thresh'))])
    colnames(thresh_df)=c('Coverage Bin', 'Habitat', 'BF Threshold')
    thresh_plot=ggplot(thresh_df, aes(x=`Coverage Bin`, y=`BF Threshold`, fill=Habitat))+
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_bar(width=0.8, position = "dodge", stat="identity")+
      labs(title=paste0(subsp, ": FPR=",100*fpr,"%")) +
      coord_cartesian(ylim = c(0.95*min(thresh_df$`BF Threshold`),NA)) + 
      scale_discrete_manual(name="", 
                            values = habitat_col, 
                            aesthetics = c("colour", "fill"))
      print(thresh_plot)
      # Make data frame of candidates
      betai_cand[[fpr_pct]][[subsp]]=betai[[subsp]][betai[[subsp]]$`BF(dB).median` > betai[[subsp]][,paste0(fpr_pct, '_FPR_BF_thresh')],]
      # Get total number of SNPs
      tot_n_snps=length(unique(betai[[subsp]]$chr_pos))
      # Make sure covariables are alphabetically ordered
      covs=unique(betai[[subsp]]$COVARIABLE_name)
      for(cov in covs[order(covs)]){
        # Get number of candidate SNPs
        cand_n_snps=length(unique(betai_cand[[fpr_pct]][[subsp]][betai_cand[[fpr_pct]][[subsp]]$COVARIABLE_name==cov,]$chr_pos))
        # FDR = number of false positives / total number of positives.
        ## number of false positives = FPR x total number SNPs.
        ## total number of positives = number of candidates.
        fdr=(fpr*tot_n_snps)/cand_n_snps
        ## FDR > 1 doesn't really make sense so I round it down to 1
        fdr=min(1, fdr)
        # Make data frame for plotting number of candidates 
        row=data.frame('Subspecies'=subsp, 
                       'Habitat'= cov, 
                       'FPR'=fpr,
                       'Total_SNPs'= tot_n_snps,
                       'Candidate_SNPs'=cand_n_snps,
                       'Expected_SNPs'=fpr*tot_n_snps,
                       'FDR'=fdr,
                       'BF_Threshold'=bf.thresh,
                       'True_Positives'=cand_n_snps*(1-fdr),
                       'False_Positives'=cand_n_snps*fdr)
        if(exists(paste0("fdr.df"))){
          fdr.df=rbind(fdr.df, row)
        }else{
          fdr.df=row
        }
      }
      # Plot to check candidates do not have biased coverage
      # Get xlim
      xmax=min(head(betai[[subsp]][order(-betai[[subsp]]$coverage),], nrow(betai[[subsp]])*top.cov)$coverage)
      coverage_plot=ggplot(NULL) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
        geom_density(data=betai[[subsp]], aes(x=coverage, col="All Sites", fill="All Sites"), alpha=0.3, bw=100) +
        geom_density(data=betai_cand[[fpr]][[subsp]], aes(x=coverage, col="Candidates", fill="Candidates"), alpha=0.3, bw=100) +
        geom_vline(xintercept = cover_bins[[subsp]], linetype='dashed', alpha=0.5) +
        labs(title=paste0(subsp, ", FPR=", fpr, ": Total Coverage Distribution"), 
             subtitle=paste0("Cropped X-axis to exclude the top ", 100*top.cov,"% coverages"), x="Total Coverage Across all Samples", y="Density") +
        scale_discrete_manual("", values = c('blue', 'red'), aesthetics = c("colour", "fill")) +
        xlim(0,xmax)
      print(coverage_plot)
    }
  }
  cat("- Number of Candidates\n")
  for(fpr in fprs){
    # Plot number positives and null expectation for a specific FPR
    ## Select relevant columns 
    fdr.bar.df=fdr.df[fdr.df$FPR==fpr, c('Subspecies', 'Habitat', 'Candidate_SNPs', 'Expected_SNPs')]
    ## Rename so plot looks nicer
    fdr.bar.df[fdr.bar.df$Subspecies=='all', 'Subspecies']="All Subspecies"
    fdr.bar.df[fdr.bar.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
    fdr.bar.df[fdr.bar.df$Subspecies=='w', 'Subspecies']="Western"
    #fdr.bar.df[fdr.bar.df$Habitat=='forest', 'Habitat']="Forest"
    #fdr.bar.df[fdr.bar.df$Habitat=='mosaic_forest', 'Habitat']="Forest-mosaic"
    #fdr.bar.df[fdr.bar.df$Habitat=='mosaic_savanna', 'Habitat']="Savanna-mosaic"
    #fdr.bar.df[fdr.bar.df$Habitat=='savanna', 'Habitat']="Savanna"
    #fdr.bar.df$Habitat=factor(fdr.bar.df$Habitat, levels=c("Forest", "Forest-mosaic", "Savanna-mosaic", "Savanna"))
    #names(habitat_col)=unique(fdr.bar.df$Habitat[order(fdr.bar.df$Habitat)])
    #fdr.bar.df$Subspecies=factor(fdr.bar.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
    ## Plot
    ### Format file for input
    colnames(fdr.bar.df)=c('Subspecies', 'Habitat', 'Number of Candidate SNPs', 'Expected Number of Candidate SNPs')
    barplot=ggplot(fdr.bar.df, aes(y=`Number of Candidate SNPs`, x=Habitat)) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.text.x = element_text(face="bold"),
            strip.background = element_rect(fill="white")) +
      geom_bar(aes(fill=Habitat), stat="identity", col='black') +
      geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`, col="Null Expectation"),lty=2)+
      #geom_errorbar(aes(y = `Expected Number of Candidate SNPs`, ymin = `Expected Number of Candidate SNPs`, ymax = `Expected Number of Candidate SNPs`), col='red') + 
      facet_wrap(~Subspecies) +
      labs(title=paste0("Number of Candidate SNPs"),
           subtitle=paste0("FPR = ", 100*fpr, "%")) +
      scale_discrete_manual(name='Habitats', values = habitat_col, aesthetics = c("fill")) +
      scale_discrete_manual(name='', values = c("Null Expectation"='red'), aesthetics = c("col")) +
      scale_x_discrete(guide = guide_axis(angle = 45)) 
    print(barplot)
  }
  return(betai_cand)
}

select_candidates.v2=function(betai, null, fprs=c(0.005, 0.0005, 0.00005), n_bins=15, top.cov=0.025, habitat_col, tail_bin_size=0.05){
  betai_cand=list()
  cover_bins=list()
  subsps=names(betai)
  for(subsp in subsps){
    cat("-", subsp, "\n")
    # Assign SNPs to coverage bins
    # Get all coverage values
    cover=betai[[subsp]]$coverage
    # Ensure that the lowest and highest bins contain 5% of the data points each
    low=max(head(cover[order(cover)], tail_bin_size*length(cover)))
    high=min(head(cover[order(-cover)], tail_bin_size*length(cover)))
    ## Separate into n bins between the low and high bin thresholds so we end up with n+2 bins over all (when you include the highest and lowest bins)
    cover_bins[[subsp]]=c(0, seq(from = low, to = high, by = (high-low)/(n_bins-2)), max(cover))
    # Assign each SNP to a bin
    betai[[subsp]]$coverage_bin=NA
    for(i in 1:(length(cover_bins[[subsp]])-1)){
      betai[[subsp]][betai[[subsp]]$coverage > cover_bins[[subsp]][i] & betai[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
      null[[subsp]][null[[subsp]]$coverage > cover_bins[[subsp]][i] & null[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
    }
    # Get BF Threshold per Bin
    for(fpr in fprs){
      fpr_pct=paste0(100*fpr, "pct")
      # Make chr_pos column
      betai[[subsp]]$chr_pos=paste(betai[[subsp]]$chr, betai[[subsp]]$pos, sep="_")
      null[[subsp]]$chr_pos=paste(null[[subsp]]$chr, null[[subsp]]$pos, sep="_")
      # Get list of habitat names
      habs=unique(betai[[subsp]]$COVARIABLE_name)
      habs=habs[order(habs)]
      # Prime the column for BF threshold per bin per habitat type for a particular FPR
      betai[[subsp]]$tmp=0
      colnames(betai[[subsp]])[colnames(betai[[subsp]])=='tmp']=paste0(fpr_pct, '_FPR_BF_thresh')
      # Loop over each bin
      bins=unique(betai[[subsp]]$coverage_bin)
      for(bin in bins[order(bins)]){
        ### Extract SNPs from the bin in the random runs (and combine the information across random runs)
        for(hab in habs){
          #### Get null BFs by combining all BFs from SNPs in the coverage bin for the given habitat type
          null_BFs=null[[subsp]][(null[[subsp]]$coverage_bin == bin)
                                 & (null[[subsp]]$COVARIABLE_name==hab), 'BF(dB).median']
          ## Get BF threshold for given fpr
          ### If n SNPs * FPR < 1 then we have to just take the highest value from the null, 
          #### I chose to report when this happens as we are no longer technically reporting a true FPR i.e. the minimum FPR we can use is determined by the number of SNPs in the null
          #### This should only be an issue if the null has very few SNPs and/or we are looking at the far extreme tail
          top_n=length(null_BFs)*fpr
          bf.thresh=min(head(null_BFs[order(-null_BFs)], max(1, top_n)))
          ## Add thresholds to dataframe
          betai[[subsp]][betai[[subsp]]$coverage_bin==bin & betai[[subsp]]$COVARIABLE_name==hab, paste0(fpr_pct, '_FPR_BF_thresh')]=bf.thresh
          ## Make dataframe for plotting
          row=data.frame('subsp'=subsp, 'COVARIABLE_name'=hab, 'fpr'=fpr, 'coverage_bin'=bin, 
                         'min'=cover_bins[[subsp]][bin], 'max'=cover_bins[[subsp]][bin+1],
                         'BF_thresh'=bf.thresh,
                         'n_null_over_thresh'=length(null_BFs[null_BFs>bf.thresh]), 'n_exome_over_thresh'=nrow(betai[[subsp]][(betai[[subsp]]$coverage_bin == bin) & 
                                                                                                (betai[[subsp]]$COVARIABLE_name==hab) &
                                                                                                betai[[subsp]]$`BF(dB).median`>bf.thresh,]), 
                         'shaded'=ifelse(top_n<1, TRUE, FALSE))
          if(exists(paste0("bin_df"))){
            bin_df=rbind(bin_df, row)
          }else{
            bin_df=row
          }
        }
      }
      # Plot number of sites in each bin
      n_sites_plot_df=unique(bin_df[bin_df$subsp==subsp & bin_df$fpr==fpr, c('coverage_bin', 'n_null_over_thresh', 'n_exome_over_thresh')])
      n_sites_plot_df=melt(n_sites_plot_df, id.vars = c("coverage_bin"))
      n_sites_plot=ggplot(n_sites_plot_df, aes(x=coverage_bin, y=value, fill=variable))+
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_bar(width=0.8, position = "dodge", stat="identity", col='black')+
        labs(title=paste0(subsp, ": Number of SNPs above FPR=",100*fpr,"% BF threshold"),
             x="Coverage bin", y="Count") +
        scale_discrete_manual(name="", 
                              labels = c("Null", "Exome"),
                              values = c('n_null_over_thresh'='blue', 'n_exome_over_thresh'='red'), 
                              aesthetics = "fill")
      print(n_sites_plot)
      # Plot thresholds
      thresh_df=bin_df[bin_df$subsp==subsp, c('coverage_bin', 'COVARIABLE_name', 'BF_thresh')]
        #unique(betai[[subsp]][c('coverage_bin', 'COVARIABLE_name', paste0(fpr_pct, '_FPR_BF_thresh'))])
      colnames(thresh_df)=c('Coverage Bin', 'COVARIABLE_name', 'BF Threshold')
      thresh_plot=ggplot(thresh_df, aes(x=`Coverage Bin`, y=`BF Threshold`, fill=COVARIABLE_name))+
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_bar(width=0.8, position = "dodge", stat="identity")+
        labs(title=paste0(subsp, ": FPR=",100*fpr,"%")) +
        coord_cartesian(ylim = c(0.95*min(thresh_df$`BF Threshold`), NA)) + 
        geom_hline(aes(yintercept=20), col='red',lty=2)+
        scale_discrete_manual(name="",
                              values = habitat_col, 
                              aesthetics = c("colour", "fill"))
      print(thresh_plot)
      # Make data frame of candidates
      betai_cand[[fpr_pct]][[subsp]]=betai[[subsp]][betai[[subsp]]$`BF(dB).median` > betai[[subsp]][,paste0(fpr_pct, '_FPR_BF_thresh')],]
      # Get total number of SNPs
      tot_n_snps=length(unique(betai[[subsp]]$chr_pos))
      # Make sure covariables are alphabetically ordered
      covs=unique(betai[[subsp]]$COVARIABLE_name)
      for(cov in covs[order(covs)]){
        # Get number of candidate SNPs
        cand_n_snps=length(unique(betai_cand[[fpr_pct]][[subsp]][betai_cand[[fpr_pct]][[subsp]]$COVARIABLE_name==cov,]$chr_pos))
        # FDR = number of false positives / total number of positives.
        ## number of false positives = FPR x total number SNPs.
        ## total number of positives = number of candidates.
        fdr=(fpr*tot_n_snps)/cand_n_snps
        ## FDR > 1 doesn't really make sense so I round it down to 1
        fdr=min(1, fdr)
        # Make data frame for plotting number of candidates 
        row=data.frame('Subspecies'=subsp, 
                       'COVARIABLE_name'= cov, 
                       'FPR'=fpr,
                       'Total_SNPs'= tot_n_snps,
                       'Candidate_SNPs'=cand_n_snps,
                       'Expected_SNPs'=fpr*tot_n_snps,
                       'FDR'=fdr,
                       'BF_Threshold'=bf.thresh,
                       'True_Positives'=cand_n_snps*(1-fdr),
                       'False_Positives'=cand_n_snps*fdr)
        if(exists(paste0("fdr.df"))){
          fdr.df=rbind(fdr.df, row)
        }else{
          fdr.df=row
        }
      }
      # Plot to check candidates do not have biased coverage
      ## Get xlim
      xmax=min(head(betai[[subsp]][order(-betai[[subsp]]$coverage),], nrow(betai[[subsp]])*top.cov)$coverage)
      ## plot
      coverage_plot=ggplot(NULL) + xlim(0,xmax)
      ## Add shaded bins? - these refer to bins where there were not enough null data points to predict such a stringent FPR, the threshold is therefore just the maximum null BF in the bin
      bin_df.tmp=bin_df[bin_df$subsp==subsp & bin_df$COVARIABLE_name==cov & bin_df$fpr==fpr & bin_df$shaded==TRUE,]
      if(nrow(bin_df.tmp)>0){
        ## Shaded areas aren't shown if their x max is beyond the graph xmax to I bring it down for the purpose of plotting
        bin_df.tmp[bin_df.tmp$max>xmax, 'max']=xmax
        coverage_plot=coverage_plot + geom_rect(data=bin_df.tmp,
                                                aes(xmin=min, xmax=max, ymin=-Inf, ymax=Inf), fill="black", alpha=0.2) 
      }
      coverage_plot=coverage_plot + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
        geom_density(data=betai[[subsp]], aes(x=coverage, col="All Sites", fill="All Sites"), alpha=0.3, bw=100) +
        geom_density(data=betai_cand[[fpr_pct]][[subsp]], aes(x=coverage, col="Candidates", fill="Candidates"), alpha=0.3, bw=100) +
        geom_vline(xintercept = cover_bins[[subsp]], linetype='dashed', alpha=0.5) +
        labs(title=paste0(subsp, ", FPR=", fpr, ": Total Coverage Distribution"), 
             subtitle=paste0("Cropped X-axis to exclude the top ", 100*top.cov,"% coverages"), x="Total Coverage Across all Samples", y="Density")
      if(length(unique(betai[[subsp]]$chr_pos))==length(unique(null[[subsp]]$chr_pos))){
        coverage_plot=coverage_plot +
          scale_discrete_manual("", values = c('blue', 'red'), aesthetics = c("colour", "fill"))
      }else{
        coverage_plot=coverage_plot +
          geom_density(data=null[[subsp]], aes(x=coverage, col="Null", fill="Null"), alpha=0.3, bw=100) +
          scale_discrete_manual("", values = c('blue', 'red', 'green'), aesthetics = c("colour", "fill"))
      }
      print(coverage_plot)
    }
  }
  cat("- Number of Candidates\n")
  for(fpr in fprs){
    # Plot number positives and null expectation for a specific FPR
    ## Select relevant columns 
    fdr.bar.df=fdr.df[fdr.df$FPR==fpr, c('Subspecies', 'COVARIABLE_name', 'Candidate_SNPs', 'Expected_SNPs')]
    ## Rename so plot looks nicer
    fdr.bar.df[fdr.bar.df$Subspecies=='all', 'Subspecies']="All Subspecies"
    fdr.bar.df[fdr.bar.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
    fdr.bar.df[fdr.bar.df$Subspecies=='w', 'Subspecies']="Western"
    #fdr.bar.df[fdr.bar.df$Habitat=='forest', 'Habitat']="Forest"
    #fdr.bar.df[fdr.bar.df$Habitat=='mosaic_forest', 'Habitat']="Forest-mosaic"
    #fdr.bar.df[fdr.bar.df$Habitat=='mosaic_savanna', 'Habitat']="Savanna-mosaic"
    #fdr.bar.df[fdr.bar.df$Habitat=='savanna', 'Habitat']="Savanna"
    #fdr.bar.df$Habitat=factor(fdr.bar.df$Habitat, levels=c("Forest", "Forest-mosaic", "Savanna-mosaic", "Savanna"))
    #names(habitat_col)=unique(fdr.bar.df$Habitat[order(fdr.bar.df$Habitat)])
    #fdr.bar.df$Subspecies=factor(fdr.bar.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
    ## Plot
    ### Format file for input
    colnames(fdr.bar.df)=c('Subspecies', 'COVARIABLE_name', 'Number of Candidate SNPs', 'Expected Number of Candidate SNPs')
    barplot=ggplot(fdr.bar.df, aes(y=`Number of Candidate SNPs`, x=COVARIABLE_name)) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.text.x = element_text(face="bold"),
            strip.background = element_rect(fill="white")) +
      geom_bar(aes(fill=COVARIABLE_name), stat="identity", col='black') +
      geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`, col="Null Expectation"),lty=2)+
      #geom_errorbar(aes(y = `Expected Number of Candidate SNPs`, ymin = `Expected Number of Candidate SNPs`, ymax = `Expected Number of Candidate SNPs`), col='red') + 
      facet_wrap(~Subspecies) +
      labs(title=paste0("Number of Candidate SNPs"),
           subtitle=paste0("FPR = ", 100*fpr, "%")) +
      scale_discrete_manual(name='Covariable', values = habitat_col, aesthetics = c("fill")) +
      scale_discrete_manual(name='', values = c("Null Expectation"='red'), aesthetics = c("col")) +
      scale_x_discrete(guide = guide_axis(angle = 45)) 
    print(barplot)
  }
  return(betai_cand)
}

assign_fpr=function(betai, null, n_bins=5, tail_bin_size=0.1, verbose=T){
  out=list()
  cover_bins=list()
  subsps=names(betai)
  for(subsp in subsps){
    if(verbose){cat("-", subsp, "\n")}
    # Assign SNPs to coverage bins
    ## Get all coverage values
    cover=betai[[subsp]]$coverage
    ## Ensure that the lowest and highest bins contain 5% of the data points each
    low=max(head(cover[order(cover)], tail_bin_size*length(cover)))
    high=min(head(cover[order(-cover)], tail_bin_size*length(cover)))
    ## Separate into n bins between the low and high bin thresholds so we end up with n+2 bins over all (when you include the highest and lowest bins)
    cover_bins[[subsp]]=c(0, seq(from = low, to = high, by = (high-low)/(n_bins-2)), max(cover))
    ## Assign each SNP to a bin
    betai[[subsp]]$coverage_bin=NA
    null[[subsp]]$coverage_bin=NA
    for(i in 1:(length(cover_bins[[subsp]])-1)){
      betai[[subsp]][betai[[subsp]]$coverage > cover_bins[[subsp]][i] & betai[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
      null[[subsp]][null[[subsp]]$coverage > cover_bins[[subsp]][i] & null[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
    }
    # Loop over each bin and assign FDR values
    bins=unique(betai[[subsp]]$coverage_bin)
    bin_fdr_out=data.frame(chr=numeric(), pos=numeric(), `BF(dB).median`=numeric(), COVARIABLE_name=character(), fdr=numeric())
    for(bin in bins[order(bins)]){
      if(verbose){cat(" Bin ", bin, "/", n_bins, "\n")}
      ## Extract SNPs from the bin in the random runs (and combine the information across random runs)
      for(hab in unique(betai[[subsp]]$COVARIABLE_name)){
        ### Get null BFs 
        null_bin=null[[subsp]][(null[[subsp]]$coverage_bin==bin)
                               & (null[[subsp]]$COVARIABLE_name==hab), c('chr', 'pos', 'BF(dB).median', 'COVARIABLE_name')]
        null_bin$dataset='null'
        #### Get FPR values for each SNP
        null_bin$fpr=rank(-null_bin$`BF(dB).median`, ties.method="max")/nrow(null_bin)
        #### Add real data
        betai_bin=betai[[subsp]][(betai[[subsp]]$coverage_bin == bin)
                               & (betai[[subsp]]$COVARIABLE_name==hab), c('chr', 'pos', 'BF(dB).median', 'COVARIABLE_name')]
        betai_bin$dataset='real'
        betai_bin$fpr=NA
        total_bin=rbind(null_bin, betai_bin)
        #### Order dataset accordingly
        total_bin=total_bin[order(total_bin$`BF(dB).median`, total_bin$dataset),]
        #### Loop over each row
        ##### The first value should have an FPR of 1
        total_bin[1, 'fpr']=1
        for(row in 2:nrow(total_bin)){
          if(is.na(total_bin[row, 'fpr'])){
            total_bin[row, 'fpr']=total_bin[row-1, 'fpr']
          }
        }
        #### Extract real data and add to file
        bin_fdr_out=rbind(bin_fdr_out, total_bin[total_bin$dataset=='real',])
      }
    }
    # Add FDR values to dataframe    
    out[[subsp]]=merge(betai[[subsp]][,c(c('MRK','chr', 'pos', 'BF(dB).median', 'M_Beta.median', 'COVARIABLE_name', 'coverage', 'coverage_bin'))], bin_fdr_out, by=c('chr', 'pos', 'BF(dB).median', 'COVARIABLE_name'), all.x=TRUE)
    out[[subsp]]=out[[subsp]][,colnames(out[[subsp]])!='dataset']
  }
  return(out)
}

assign_fpr.v2=function(betai, null, n_bins=5, tail_bin_size=0.1, beta_split=F, verbose=T){
  out=list()
  cover_bins=list()
  subsps=names(betai)
  for(subsp in subsps){
    if(verbose){cat("-", subsp, "\n")}
    # Assign SNPs to coverage bins
    ## Get all coverage values
    cover=betai[[subsp]]$coverage
    ## Ensure that the lowest and highest bins contain 5% of the data points each
    low=max(head(cover[order(cover)], tail_bin_size*length(cover)))
    high=min(head(cover[order(-cover)], tail_bin_size*length(cover)))
    ## Separate into n bins between the low and high bin thresholds so we end up with n+2 bins over all (when you include the highest and lowest bins)
    cover_bins[[subsp]]=c(0, seq(from = low, to = high, by = (high-low)/(n_bins-2)), max(cover))
    ## Assign each SNP to a bin
    betai[[subsp]]$coverage_bin=NA
    null[[subsp]]$coverage_bin=NA
    for(i in 1:(length(cover_bins[[subsp]])-1)){
      betai[[subsp]][betai[[subsp]]$coverage > cover_bins[[subsp]][i] & betai[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
      null[[subsp]][null[[subsp]]$coverage > cover_bins[[subsp]][i] & null[[subsp]]$coverage <= cover_bins[[subsp]][i+1], 'coverage_bin']=i
    }
    # Loop over each bin and assign FDR values
    bins=unique(betai[[subsp]]$coverage_bin)
    bin_fpr_out=data.frame(chr=numeric(), pos=numeric(), `BF(dB).median`=numeric(), M_Beta.median=numeric(), COVARIABLE_name=character(), fpr=numeric())
    for(bin in bins[order(bins)]){
      ## Extract SNPs from the bin in the random runs (and combine the information across random runs)
      for(hab in unique(betai[[subsp]]$COVARIABLE_name)){
        ### Get null BFs 
        null_bin=null[[subsp]][(null[[subsp]]$coverage_bin==bin)
                               & (null[[subsp]]$COVARIABLE_name==hab), c('chr', 'pos', 'BF(dB).median', 'M_Beta.median','COVARIABLE_name')]
        null_bin$dataset='null'
        #### Get FPR values for each SNP
        null_bin$fpr=rank(-null_bin$`BF(dB).median`, ties.method="max")/nrow(null_bin)
        #### Add real data
        betai_bin=betai[[subsp]][(betai[[subsp]]$coverage_bin == bin)
                                 & (betai[[subsp]]$COVARIABLE_name==hab), c('chr', 'pos', 'BF(dB).median', 'M_Beta.median','COVARIABLE_name')]
        betai_bin$dataset='real'
        betai_bin$fpr=NA
        total_bin=rbind(null_bin, betai_bin)
      }
      # Estimate FPR independently for sites with positive and negative correlation coeffeincts to account for biases???
      if(beta_split){
        total_bin.l=list()
        total_bin.l[['pos']]=total_bin[total_bin$M_Beta.median>0,]
        total_bin.l[['neg']]=total_bin[total_bin$M_Beta.median<=0,]
        for(dir in names(total_bin.l)){
          #### Order dataset accordingly
          total_bin.l[[dir]]=total_bin.l[[dir]][order(total_bin.l[[dir]]$`BF(dB).median`, total_bin.l[[dir]]$dataset),]
          #### Loop over each row
          ##### The first value should have an FPR of 1
          total_bin.l[[dir]][1, 'fpr']=1
          for(row in 2:nrow(total_bin.l[[dir]])){
            if(is.na(total_bin.l[[dir]][row, 'fpr'])){
              total_bin.l[[dir]][row, 'fpr']=total_bin.l[[dir]][row-1, 'fpr']
            }
          }
          total_bin=rbind(total_bin.l[['pos']], total_bin.l[['neg']])
        }
      }else{
        #### Order dataset accordingly
        total_bin=total_bin[order(total_bin$`BF(dB).median`, total_bin$dataset),]
        #### Loop over each row
        ##### The first value should have an FPR of 1
        total_bin[1, 'fpr']=1
        for(row in 2:nrow(total_bin)){
          if(is.na(total_bin[row, 'fpr'])){
            total_bin[row, 'fpr']=total_bin[row-1, 'fpr']
          }
        }
      }
      #### Extract real data and add to file
      bin_fpr_out=rbind(bin_fpr_out, total_bin[total_bin$dataset=='real',])
    }
    # Add FPR values to data frame    
    out[[subsp]]=merge(betai[[subsp]][,c('MRK','chr', 'pos', 'BF(dB).median', 'M_Beta.median', 'COVARIABLE_name', 'coverage', 'coverage_bin')], bin_fpr_out, by=c('chr', 'pos', 'BF(dB).median', 'M_Beta.median', 'COVARIABLE_name'), all.x=TRUE)
    out[[subsp]]=out[[subsp]][,colnames(out[[subsp]])!='dataset']
  }
  return(out)
}

plot_fpr_stats=function(betai, null=betai_chr21, fprs=c(0.005, 0.001, 0.0005), top.cov=0.025, N_SNPs=TRUE, BF_thresh=TRUE, coverage=TRUE){
  subsp=names(betai)
  fdr.df=NULL
  cat("Coverage Bin Stats\n")
  for(subsp in subsps){
    if(subsp=='all'){
        subsp_name="All Subspecies"
      }else if(subsp=='ce'){
        subsp_name="Central-Eastern"
      }else if(subsp=='w'){
        subsp_name="Western"
      }else{
        subsp_name=subsp
        }
    for(fpr in c(fprs)){
      bin_df=data.frame(fpr=numeric(), bin=numeric(), COVARIABLE_name=character(), n_SNPs=numeric(), dataset=character(), bf.thresh=numeric())
      # Get total number of SNPs
      tot_n_snps=nrow(unique(betai[[subsp]][,c('chr', 'pos')]))
      # Make sure covariables are alphabetically ordered
      covs=unique(betai[[subsp]]$COVARIABLE_name)
      # For each covariable
      for(cov in covs[order(covs)]){
        # Select candidates
        cands=unique(betai[[subsp]][betai[[subsp]]$fpr<=fpr & betai[[subsp]]$COVARIABLE_name==cov, ])
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
                       'COVARIABLE_name'= cov, 
                       'FPR'=fpr,
                       'Total_SNPs'= tot_n_snps,
                       'Candidate_SNPs'=cand_n_snps,
                       'Expected_SNPs'=fpr*tot_n_snps,
                       'BF_Threshold'= min(cands$`BF(dB).median`),
                       'FDR'=fdr,
                       'True_Positives'=cand_n_snps*(1-fdr),
                       'False_Positives'=cand_n_snps*fdr)
        if(!is.null(fdr.df)){
          fdr.df=rbind(fdr.df, row)
        }else{
          fdr.df=row
        }
        # For each coverage bin
        bins=unique(betai[[subsp]]$coverage_bin)
        for(bin in bins[order(bins)]){
          # Here I calculate the number of null SNPs which were used to calculate the FPR given the FPRs
          ## FPR is rank/total SNPs in bin - we want the maximum rank (which tells us how many null SNPs in that bin are over the threshold)
          ## First calculate total number of SNPs in bin by inputting rank 1 and and the minimum fpr (which must correspond to rank 1) 
          n_null_snps_bin=1/min(cands[cands$coverage_bin==bin,'fpr'])
          ## Now we have the total number of SNPs and the maximum FPR (corresponding to the SNP with the highest rank), we can calculate the highest rank
          ### Because in tied instances, rank takes the largest number, this tells us the total number of null SNPs above this threshold
          n_null_snps_bin_over_thresh=n_null_snps_bin*max(cands[cands$coverage_bin==bin,'fpr'])
          # Get the BF threshold for the bin
          bf.thresh=min(cands[cands$coverage_bin==bin, 'BF(dB).median'])
          # Add row for null data
          bin_df=rbind(bin_df, data.frame(fpr=fpr,
                                          bin=bin, 
                                          COVARIABLE_name=cov, 
                                          n_SNPs=n_null_snps_bin_over_thresh, 
                                          dataset='Null', 
                                          bf.thresh=bf.thresh))
          # Add row for real data (having it in long format makes it easier to plot)
          bin_df=rbind(bin_df, data.frame(fpr=fpr,
                                          bin=bin, 
                                          COVARIABLE_name=cov, 
                                          n_SNPs=nrow(cands[cands$coverage_bin==bin,]), 
                                          dataset='Exome', 
                                          bf.thresh=bf.thresh))
        }
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
        thresh_plot=ggplot(bin_df[bin_df$bf.thresh!=Inf,], aes(x=bin, y=bf.thresh, fill=COVARIABLE_name))+
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5),
                panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
          geom_bar(width=0.8, position = "dodge", stat="identity", col="black")+
          labs(title=paste0(subsp_name, ", FPR<",100*fpr,"%"),
               y="BF(dB) threshold", x="Coverage bin") +
          #coord_cartesian(ylim = c(min(19.5, 0.95*min(bin_df$bf.thresh)), NA)) + 
          coord_cartesian(ylim = c(0.95*min(bin_df$bf.thresh), NA)) + 
          #geom_hline(aes(yintercept=20), col='red',lty=2)+
          scale_discrete_manual(name="",
                                values = habitat_col, 
                                aesthetics = c("colour", "fill"))
        if(length(covs)==1){
          thresh_plot=thresh_plot+guides(fill = F)
        }
        print(thresh_plot)
      }
      # Plot to check candidates do not have biased coverage
      ## Get (approximate) limits of bins using the maximum coverages observed in each bin
      cover_bins=aggregate(coverage ~ coverage_bin, betai[[subsp]], function(x) max(x))$coverage
      ## Get xmax (to ignore the top x% of coverage values and make the graph readable)
      xmax=min(head(betai[[subsp]][order(-betai[[subsp]]$coverage),], nrow(betai[[subsp]])*top.cov)$coverage)
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
                legend.background = element_rect(fill='transparent')) +
          #geom_density(data=betai[[subsp]], aes(x=coverage, col="Exome", fill="Exome"), alpha=0.3, bw=100) +
          #geom_density(data=cands, aes(x=coverage, col="Candidates", fill="Candidates"), alpha=0.3, bw=100) +
          #geom_density(data=null[[subsp]], aes(x=coverage, col="Non-genic\n-chr21", fill="Non-genic\n-chr21"), alpha=0.3, bw=100) +
          stat_bin(data=betai[[subsp]], aes(x=coverage, y=..density.., col="Exome", fill="Exome"), alpha=1, binwidth=binwidth, geom="step", size=lw,
                   position=position_nudge(x=-0.5*binwidth)) +
          geom_histogram(data=betai[[subsp]], aes(x=coverage, y=..density.., fill="Exome"), alpha=alpha, binwidth=binwidth, size=0) +
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
          #scale_discrete_manual("", values = c('black', 'turquoise4', 'violetred'), aesthetics = c("colour", "fill")) + 
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
    fdr.bar.df=fdr.df[fdr.df$FPR==fpr, c('Subspecies', 'COVARIABLE_name', 'Candidate_SNPs', 'Expected_SNPs')]
    if(nrow(fdr.bar.df)>0){
      ## Rename so plot looks nicer
      fdr.bar.df[fdr.bar.df$Subspecies=='all', 'Subspecies']="All Subspecies"
      fdr.bar.df[fdr.bar.df$Subspecies=='ce', 'Subspecies']="Central-Eastern"
      fdr.bar.df[fdr.bar.df$Subspecies=='w', 'Subspecies']="Western"
      #fdr.bar.df$Subspecies=factor(fdr.bar.df$Subspecies, levels=c("All Subspecies", "Central-Eastern", "Western"))
      ## Plot
      ### Format file for input
      colnames(fdr.bar.df)=c('Subspecies', 'COVARIABLE_name', 'Number of Candidate SNPs', 'Expected Number of Candidate SNPs')
      barplot=ggplot(fdr.bar.df, aes(y=`Number of Candidate SNPs`, x=COVARIABLE_name)) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              strip.text.x = element_text(face="bold"),
              strip.background = element_rect(fill="white"),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        geom_bar(aes(fill=COVARIABLE_name), stat="identity", col='black') +
        #geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`, col="Null Expectation"),lty=2, size=1)+
        geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`), col="white", size=3)+
        geom_hline(aes(yintercept=`Expected Number of Candidate SNPs`, col="Null Expectation"), size=2)+
        #geom_errorbar(aes(y = `Expected Number of Candidate SNPs`, ymin = `Expected Number of Candidate SNPs`, ymax = `Expected Number of Candidate SNPs`), col='red') + 
        facet_wrap(~Subspecies) +
        labs(title=paste0("Number of Candidate SNPs"),
             subtitle=paste0("FPR = ", 100*fpr, "%")) +
        scale_discrete_manual(name='Covariable', values = habitat_col, aesthetics = c("fill")) +
        scale_discrete_manual(name='', values = c("Null Expectation"='violetred'), aesthetics = c("col")) +
        scale_x_discrete(guide = guide_axis(angle = 45)) 
      print(barplot)
    }
  }
}

candidate_overlap=function(betai_cand){
  # Overlap between categories
  for(fpr in names(betai_cand)){
    for(subsp in names(betai_cand[[fpr]])){
      venn_list=list()
      covs=unique(betai_cand[[fpr]][[subsp]]$COVARIABLE_name)
      for(cov in covs){
        venn_list[[cov]]=betai_cand[[fpr]][[subsp]][betai_cand[[fpr]][[subsp]]$COVARIABLE_name==cov, 'chr_pos']
      }
      if(length(venn_list)>1){
        venn=ggVennDiagram(venn_list) +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
          labs(title=paste0(subsp), subtitle=paste0("FPR=", fpr)) +
          scale_fill_gradient(low="blue", high = "red")
        print(venn)
      }
    }
  }
  # Overlap between subspecies datasets
  if(length(subsps>1)){
    for(fpr in names(betai_cand)){
      covs=unique(betai_cand[[fpr]][[1]]$COVARIABLE_name)
      for(cov in covs){
        venn_list=list()
        for(subsp in names(betai_cand[[fpr]])){
          venn_list[[subsp]]=betai_cand[[fpr]][[subsp]][betai_cand[[fpr]][[subsp]]$COVARIABLE_name==cov, 'chr_pos']
        }
        if(length(venn_list)>1){
          venn=ggVennDiagram(venn_list) +
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
            labs(title=paste0(cov), subtitle=paste0("FPR=", fpr)) +
            scale_fill_gradient(low="blue",high = "red")
          print(venn)
        }
      }
    }
  }
}

candidate_overlap.v2=function(betai, fprs=c(0.005, 0.001, 0.0005)){
  subsps=names(betai)
  # Overlap between categories
  for(fpr in fprs){
    for(subsp in subsps){
      venn_list=list()
      covs=unique(betai[[subsp]]$COVARIABLE_name)
      for(cov in covs){
        snps=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, c('chr', 'pos')]
        venn_list[[cov]]=paste(snps$chr, snps$pos, sep="_")
      }
      if(length(venn_list)>1){
        venn=ggVennDiagram(venn_list) +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
          labs(title=paste0(subsp), subtitle=paste0("FPR=", fpr)) +
          scale_fill_gradient(low="blue", high = "red")
        print(venn)
      }
    }
  }
  # Overlap between subspecies datasets
  if(length(subsps>1)){
    for(fpr in fprs){
      covs=unique(betai[[1]]$COVARIABLE_name)
      for(cov in covs){
        venn_list=list()
        for(subsp in subsps){
          snps=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, c('chr', 'pos')]
          venn_list[[subsp]]=paste(snps$chr, snps$pos, sep="_")
        }
        if(length(venn_list)>1){
          venn=ggVennDiagram(venn_list) +
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
            labs(title=paste0(cov), subtitle=paste0("FPR=", fpr)) +
            scale_fill_gradient(low="blue",high = "red")
          print(venn)
        }
      }
    }
  }
}

candidate_overlap.v3=function(betai, fprs=c(0.005, 0.001, 0.0005)){
  subsps=names(betai)
  # Overlap between categories
  for(fpr in fprs){
    for(subsp in subsps){
      venn_list=list()
      covs=unique(betai[[subsp]]$COVARIABLE_name)
      for(cov in covs){
        snps=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, c('chr', 'pos')]
        venn_list[[cov]]=paste(snps$chr, snps$pos, sep="_")
      }
      if(length(venn_list)>1){
        venn=ggVennDiagram(venn_list) +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
          labs(subtitle=paste0(subsp), title=paste0("Candidate SNPs: FPR < ", fpr)) +
          scale_fill_viridis()
        print(venn)
      }
    }
  }
  ## If the file has gene annotations, look at candidate gene overlap
  if('gene' %in% colnames(betai[[1]])){
    for(fpr in fprs){
      for(subsp in subsps){
        venn_list=list()
        covs=unique(betai[[subsp]]$COVARIABLE_name)
        for(cov in covs){
          genes=unique(betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, 'gene'])
          venn_list[[cov]]=genes
        }
        if(length(venn_list)>1){
          venn=ggVennDiagram(venn_list) +
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
            labs(subtitle=paste0(subsp), title=paste0("Candidate Genes: FPR < ", fpr))  +
            scale_fill_viridis()
          print(venn)
        }
      }
    }
  }
  # Overlap between subspecies datasets
  if(length(subsps>1)){
    for(fpr in fprs){
      covs=unique(betai[[1]]$COVARIABLE_name)
      for(cov in covs){
        venn_list=list()
        for(subsp in subsps){
          snps=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, c('chr', 'pos')]
          venn_list[[subsp]]=paste(snps$chr, snps$pos, sep="_")
        }
        if(length(venn_list)>1){
          venn=ggVennDiagram(venn_list) +
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
            labs(subtitle=paste0(cov), title=paste0("Candidate SNPs: FPR < ", fpr)) +
            scale_fill_viridis()
          print(venn)
        }
      }
    }
    ## If the file has gene annotations, look at candidate gene overlap
    if('gene' %in% colnames(betai[[1]])){
      for(fpr in fprs){
        covs=unique(betai[[1]]$COVARIABLE_name)
        for(cov in covs){
          venn_list=list()
          for(subsp in subsps){
            genes=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, 'gene']
            venn_list[[subsp]]=genes
          }
          if(length(venn_list)>1){
            venn=ggVennDiagram(venn_list) +
              theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
              labs(subtitle=paste0(cov), title=paste0("Candidate Genes: FPR < ", fpr))  +
              scale_fill_viridis()
            print(venn)
          }
        }
      }
    }
  }
}

candidate_overlap.v4=function(betai, fprs=c(0.005, 0.001, 0.0005), set_color=c("All Samples"="black", "Central-Eastern"="brown4", "Nigeria-Cameroon"="red", "Western"="blue"), title_suffix=""){
  subsps=names(betai)
  # Venn function
  plot_venn=function(venn_list, title="", subtitle="", set_color=set_color){
    plot=ggVennDiagram(venn_list) +
      theme(plot.title = element_text(hjust = 0.5, size=15), plot.subtitle = element_text(hjust = 0.5, size=10)) +
      labs(subtitle=subtitle, title=title) +
      scale_fill_viridis() +
      #scale_color_manual(values = set_color) +
      scale_x_continuous(expand = expansion(mult = 0.2))
    print(plot)
  }
  # Overlap between categories
  for(fpr in fprs){
    for(subsp in subsps){
      venn_list=list()
      covs=unique(betai[[subsp]]$COVARIABLE_name)
      for(cov in covs){
        snps=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, c('chr', 'pos')]
        venn_list[[cov]]=paste(snps$chr, snps$pos, sep="_")
      }
      if(length(venn_list)>1){
        plot_venn(venn_list, title=paste0("Candidate SNPs: FPR < ", 100*fpr, "%", title_suffix), subtitle=subsp)
      }
    }
  }
  ## If the file has gene annotations, look at candidate gene overlap
  if('gene' %in% colnames(betai[[1]])){
    for(fpr in fprs){
      for(subsp in subsps){
        venn_list=list()
        covs=unique(betai[[subsp]]$COVARIABLE_name)
        for(cov in covs){
          genes=unique(betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, 'gene'])
          venn_list[[cov]]=genes
        }
        if(length(venn_list)>1){
          plot_venn(venn_list, title=paste0("Candidate Genes: FPR < ", 100*fpr, "%", title_suffix), subtitle=subsp)
        }
      }
    }
  }
  # Overlap between subspecies datasets
  if(length(subsps>1)){
    for(fpr in fprs){
      covs=unique(betai[[1]]$COVARIABLE_name)
      for(cov in covs){
        venn_list=list()
        for(subsp in subsps){
          snps=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, c('chr', 'pos')]
          venn_list[[subsp]]=paste(snps$chr, snps$pos, sep="_")
        }
        if(length(venn_list)>1){
          plot_venn(venn_list, title=paste0("Candidate SNPs: FPR < ", 100*fpr, "%", title_suffix), subtitle=cov)
        }
      }
    }
    ## If the file has gene annotations, look at candidate gene overlap
    if('gene' %in% colnames(betai[[1]])){
      for(fpr in fprs){
        covs=unique(betai[[1]]$COVARIABLE_name)
        for(cov in covs){
          venn_list=list()
          for(subsp in subsps){
            genes=betai[[subsp]][betai[[subsp]]$COVARIABLE_name==cov & betai[[subsp]]$fpr<=fpr, 'gene']
            venn_list[[subsp]]=genes
          }
          if(length(venn_list)>1){
            plot_venn(venn_list, title=paste0("Candidate Genes: FPR < ", 100*fpr, "%", title_suffix), subtitle=cov)
          }
        }
      }
    }
  }
}

snps_per_gene_across_tails=function(tails_list, select_tails=NULL){
  subsps=names(tails_list)
  for(subsp in subsps){
    cat(subsp, "\n")
    # Create empty output data frame
    out.df=data.frame(tail=numeric(), SNPs.per.genes_data=numeric(), SNPs.per.genes_null=numeric())
    # Prepare background snps file
    bg_snps=unique(betai.ann[[subsp]][,c('chr', 'pos', 'gene'),])
    bg_snps$chr_pos=paste(bg_snps$chr, bg_snps$pos, sep="_")
    # Select specific tails?
    if(!is.null(select_tails)){
      tails=select_tails
    }else{
      tails=names(tails_list[[subsp]])
    }
    for(tail in tails){
      tails_list[[subsp]][[tail]]$chr_pos=paste(tails_list[[subsp]][[tail]]$chr, tails_list[[subsp]][[tail]]$pos, sep="_")
      for(cov in unique(tails_list[[subsp]][[tail]]$COVARIABLE_name)){
        cands_tail=tails_list[[subsp]][[tail]][tails_list[[subsp]][[tail]]$COVARIABLE_name==cov,]
        # Only include SNPs in genes
        cands_tail=cands_tail[!is.na(cands_tail$gene),]
        if(nrow(cands_tail)>0){
          # Prepare null df
          null_df=data.frame(chr=numeric(), pos=numeric(), gene=character())
          # Ensure we have unique rows
          cands_tail=unique(cands_tail)
          ### Get the number of tail SNPs
          n_snps=length(unique(cands_tail$chr_pos))
          ### Get the number of genes which overlap tail SNPs
          n_genes=length(unique(cands_tail$gene))
          ## Null
          n_genes_null=c()
          ### Do the following 10 times and use the mean
          for(i in 1:10){
            #### Randomly sample the same number of SNPs in the tail from all sites
            rand.snps=unique(bg_snps$chr_pos)[sample(length(unique(bg_snps$chr_pos)), n_snps)]
            rand.snps=bg_snps[bg_snps$chr_pos %in% rand.snps,]
            #### Get number of genes which overlap with randomly sampled SNPs
            n_genes_null[i]=length(unique(rand.snps$gene))
            #### Add suffix to gene names (so results from each random run are unique)
            rand.snps$gene=paste0(rand.snps$gene, paste0(" ",i))
            null_df=rbind(null_df, rand.snps)
          }
          #### Get the mean
          n_genes_null=mean(n_genes_null)
          ## Add to output
          out.df=rbind(out.df, data.frame(tail=tail, SNPs.per.genes_data=n_snps/n_genes, SNPs.per.genes_null=n_snps/n_genes_null))
          # Plot genes per snp
          ## Real data
          snps_per_gene.data=data.frame(table(cands_tail$gene))
          colnames(snps_per_gene.data)=c('SNPs per Gene', 'Frequency')
          snps_per_gene.data$Dataset='Data'
          ### get xmax
          ### xmax is set as the maximum value obsevred in the real data (as we take more samples to geenrate the null we get a small humber of high values which messes uop the scale)
          xmax=max(as.numeric(as.character(data.frame(table(snps_per_gene.data$Frequency))$Var1)))+0.5
          ## Null
          snps_per_gene.null=data.frame(table(null_df$gene))
          colnames(snps_per_gene.null)=c('SNPs per Gene', 'Frequency')
          snps_per_gene.null$Dataset='Null'
          ## Combine
          snps_per_gene=rbind(snps_per_gene.data, snps_per_gene.null)
          hist=ggplot(snps_per_gene, aes(x=Frequency, fill=Dataset)) +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
            geom_histogram(aes(y = ..density..), position = 'identity', alpha=0.5, bins=xmax+1) +
            scale_discrete_manual(name='', values = c('Data' = 'red', 'Null' = 'blue'), aesthetics='fill') +
            labs(title="SNPs per Gene Histogram",
                 subtitle=paste0("Subspecies dataset: ", subsp, ", tail: ",tail,", Covariable: " , cov),
                 x="SNPs per Gene",
                 y = "Density") +
            xlim(c(0.5, xmax))
          print(hist)
        }
      }
    }
    # Ensure x axis is ordered according to the order of the input list
    out.df$tail=factor(out.df$tail, levels=names(tails_list[[subsp]]))
    # Plot
    plot=ggplot(out.df) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      geom_point(aes(x=tail, y=SNPs.per.genes_data, col="Data")) +
      geom_point(aes(x=tail, y=SNPs.per.genes_null, col="Null")) +
      scale_discrete_manual(name='', values = c('Data' = 'red', 'Null' = 'blue'), aesthetics='colour') +
      labs(title="Genes per SNP across tails\ncompared to random null",
           subtitle=paste0("Subspecies dataset: ", subsp, ", Covariable: " , cov),
           x="Tail",
           y = "# SNPs/# Genes")
    print(plot)
  }
}

snps_per_gene_across_tails.v2=function(betai_fpr, fprs, reps=10){
  # Create empty output data frame
  out.df=data.frame(tail=numeric(), SNPs.per.genes_data=numeric(), SNPs.per.genes_null=numeric())
  # Prep unique snp identifier column (chr_pos)
  betai_fpr$chr_pos=paste(betai_fpr$chr, betai_fpr$pos, sep="_")
  # Prepare background snps file
  bg_snps=unique(betai_fpr[,c('chr', 'pos', 'gene'),])
  bg_snps$chr_pos=paste(bg_snps$chr, bg_snps$pos, sep="_")
  # Output lists
  hist_spg=list()
  hist_gps=list()
  for(fpr in fprs){
    cands=betai_fpr[betai_fpr$fpr<=fpr,]
    # Only include SNPs in genes
    cands=cands[!is.na(cands$gene),]
    if(nrow(cands)>0){
      # Prepare null df
      null_df=data.frame(chr=numeric(), pos=numeric(), gene=character())
      # Ensure we have unique rows
      cands=unique(cands)
      ### Get the number of tail SNPs
      n_snps=length(unique(cands$chr_pos))
      ### Get the number of genes which overlap tail SNPs
      n_genes=length(unique(cands$gene))
      ## Null
      n_genes_null=c()
      ### Do the following 10 times and use the mean
      for(i in 1:10){
        #### Randomly sample the same number of SNPs in the tail from all sites
        rand.snps=unique(bg_snps$chr_pos)[sample(length(unique(bg_snps$chr_pos)), n_snps)]
        rand.snps=bg_snps[bg_snps$chr_pos %in% rand.snps,]
        #### Get number of genes which overlap with randomly sampled SNPs
        n_genes_null[i]=length(unique(rand.snps$gene))
        #### Add suffix to gene names (so results from each random run are unique)
        rand.snps$gene=paste0(rand.snps$gene, paste0(" ",i))
        null_df=rbind(null_df, rand.snps)
      }
      #### Get the mean
      n_genes_null=mean(n_genes_null)
      ## Add to output
      out.df=rbind(out.df, data.frame(fpr=fpr, SNPs.per.genes_data=n_snps/n_genes, SNPs.per.genes_null=n_snps/n_genes_null))
      
      # SNPs per gene
      ## Real data
      snps_per_gene.data=data.frame(table(cands$gene))
      colnames(snps_per_gene.data)=c('SNPs per Gene', 'Frequency')
      snps_per_gene.data$Dataset='Data'
      ## Null
      snps_per_gene.null=data.frame(table(null_df$gene))
      colnames(snps_per_gene.null)=c('SNPs per Gene', 'Frequency')
      snps_per_gene.null$Dataset='Null'
      ## Combine
      snps_per_gene=rbind(snps_per_gene.data, snps_per_gene.null)
      ### get xmax
      xmax=max(as.numeric(as.character(data.frame(table(snps_per_gene$Frequency))$Var1)))+0.5
      hist_spg[[paste(100*fpr,"pct")]]=ggplot(snps_per_gene, aes(x=Frequency, fill=Dataset)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
        geom_histogram(aes(y = ..density..), position = 'identity', alpha=0.5, bins=xmax+1) +
        scale_discrete_manual(name='', values = c('Data' = 'red', 'Null' = 'blue'), aesthetics='fill') +
        labs(title="SNPs per Gene Histogram",
             subtitle=paste0("Subspecies dataset: ", subsp, ", fpr: ",fpr,", Covariable: " , cov),
             x="SNPs per Gene",
             y = "Density") +
        xlim(c(0.5, xmax))
      
      # Genes per SNP
      ## Real data
      genes_per_snp.data=data.frame(table(cands$chr_pos))
      colnames(genes_per_snp.data)=c('Genes per SNP', 'Frequency')
      genes_per_snp.data$Dataset='Data'
      ## Null
      genes_per_snp.null=data.frame(table(null_df$chr_pos))
      colnames(genes_per_snp.null)=c('Genes per SNP', 'Frequency')
      genes_per_snp.null$Dataset='Null'
      ## Combine
      genes_per_snp=rbind(genes_per_snp.data, genes_per_snp.null)
      ### get xmax
      xmax=max(as.numeric(as.character(data.frame(table(genes_per_snp$Frequency))$Var1)))+0.5
      hist_gps[[paste(100*fpr,"pct")]]=ggplot(genes_per_snp, aes(x=Frequency, fill=Dataset)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
        geom_histogram(aes(y = ..density..), position = 'identity', alpha=0.5, bins=xmax+1) +
        scale_discrete_manual(name='', values = c('Data' = 'red', 'Null' = 'blue'), aesthetics='fill') +
        labs(title="Genes per SNP Histogram",
             subtitle=paste0("Subspecies dataset: ", subsp, ", fpr: ",fpr,", Covariable: " , cov),
             x="Genes per SNP",
             y = "Density") +
        xlim(c(0.5, xmax))
    }
  }
  for(plot in hist_spg){print(plot)}
  for(plot in hist_gps){print(plot)}
  # Ensure x axis is ordered according to the order of the input list
  out.df$fpr=factor(out.df$fpr, levels=fprs)
  # Plot
  plot=ggplot(out.df) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    geom_point(aes(x=fpr, y=SNPs.per.genes_data, col="Data")) +
    geom_point(aes(x=fpr, y=SNPs.per.genes_null, col="Null")) +
    scale_discrete_manual(name='', values = c('Data' = 'red', 'Null' = 'blue'), aesthetics='colour') +
    labs(title="SNPs per genes across tails\ncompared to random null",
         subtitle=paste0("Subspecies dataset: ", subsp, ", Covariable: " , cov),
         x="FPR",
         y = "# SNPs/# Genes")
  print(plot)
}

