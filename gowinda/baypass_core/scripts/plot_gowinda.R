plot_gowinda=function(gowinda_out, FDR.thresh=0.05, return.significant=FALSE, title="Gowinda Results"){
  gowinda_out=as.data.frame(gowinda_out)
  colnames(gowinda_out)=c("gene_set", "expected_hits", "observed_hits", "p", "FDR", "No._candidate_genes_in_gene_set", "No._genes_in_gene_set_and_annotation_and_snp", "No._genes_in_gene_set", "gene_set_description", "candidate_genes")
  # Add expected and observed column for plotting
  gowinda_out$expected_and_observed=paste0("E:",round(gowinda_out$expected_hits,2) , " O:", gowinda_out$observed_hits)
  # Plot
  ## Select top n most enriched categories
  gowinda_out_plot=head(gowinda_out[order(gowinda_out$FDR), ], 30)
  ## Trim the description so it isn't ridiculously long 
  gowinda_out_plot$plot.x=strtrim(gowinda_out_plot$gene_set_description, 30)
  ### Add '...' to the end of trimmed names
  gowinda_out_plot$plot.x=gsub('^(.{40})(.*)$', '\\1...\\2', gowinda_out_plot$plot.x)
  ### Add gene set code to ensure all x values are unique
  gowinda_out_plot$plot.x=paste0(gowinda_out_plot$gene_set, "; ",gowinda_out_plot$plot.x)
  ## Ensure bars will be in the correct order 
  gowinda_out_plot=gowinda_out_plot[order(-gowinda_out_plot$FDR),]
  gowinda_out_plot$plot.x=factor(gowinda_out_plot$plot.x, levels = gowinda_out_plot$plot.x)
  ## Make FDR rounded column for writing on plot
  gowinda_out_plot$FDR_rounded=round(gowinda_out_plot$FDR, 2)
  ## ggplot
  plot=ggplot(gowinda_out_plot, aes(x=plot.x, y=FDR, fill=FDR)) + 
    ## General
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size=9)) +
    labs(title=paste0(title), x="Gene set", y="FDR") +
    coord_flip() +
    ## Axes
    ylim(0,1.25) +
    ## Plot
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low='red', mid='yellow', high='blue', midpoint = 0.5, limit = c(0, 1)) +
    # Text
    geom_text(aes(label=FDR_rounded), hjust=1.25, color="black", size=3) +
    geom_text(aes(y=1.15,label=expected_and_observed), color="black", size=3) +
    ## FDR threshold line
    geom_hline(yintercept = FDR.thresh, linetype="dashed", colour="black")
  print(plot)
  ## Select significantly enriched gene sets
  if(return.significant==TRUE){
    return(gowinda_out[gowinda_out$FDR<=FDR.thresh,])
  }
}

read_and_plot_gowinda=function(subsps, tails, gene_set, FDR.thresh=0.05, return.significant=FALSE){
  significant=list()
  for(subsp in subsps){
    for(tail in tails){
      # Using read.table() didn't work well because for some reason it didn't recognise some tabs as deliminators 
      gowinda_out=data.frame(fread(paste0("/Users/harrisonostridge/OneDrive\ -\ University\ College\ London/Projects/PanAf/phase1and2_exomes/gowinda/baypass_core/output/gowinda_output/",subsp,"/",subsp,"_enrichment_",tail,".tail_",gene_set,".txt"), sep='\t', header=F))
      significant[[paste0(subsp,"_",tail)]]=plot_gowinda(gowinda_out, FDR.thresh, return.significant, title = paste0(subsp, " ", tail, " ", gene_set))
    }
  }
  if(return.significant==TRUE){return(significant)}
}
