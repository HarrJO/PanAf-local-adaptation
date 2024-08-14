# allele_frequency_tools.R
## General tools for analysing the output from ANGSD
# Library
library(ggplot2)

exome.vs.chr21.sfs=function(exome.ac, chr21.ac, fixed.sites=F, title="Exome and chr21 SFS", bins=100, log=TRUE){
  # Get the subspecies prefixes
  chr21_pop_names=names(chr21.ac)
  chr21_pop_names=chr21_pop_names[3:length(chr21_pop_names)]
  chr21_subsps=unique(substr(chr21_pop_names, start = 1 , stop = 1))
  exome_pop_names=names(exome.ac)
  exome_pop_names=exome_pop_names[3:length(exome_pop_names)]
  exome_subsps=unique(substr(exome_pop_names, start = 1 , stop = 1))
  subsps=unique(c(chr21_subsps, exome_subsps))
  if(length(subsps)>1){subsps=c("", subsps)}else{subsps=""}
  for(subsp in subsps){
    if(subsp=='c'){
      subtitle="Central"
    }else if(subsp=='e'){
      subtitle="Eastern"
    }else if(subsp=='n'){
      subtitle="Nigeria-Cameroon"
    }else if(subsp=='w'){
      subtitle="Western"
    }else{subtitle=""}
    if(sum(grepl(paste0("^", subsp), names(chr21.ac)))>1 & sum(grepl(paste0("^", subsp), names(exome.ac)))>1){
      # Select columns corresponding to the subspecies
      subsp.chr21.ac=chr21.ac[, grepl(paste0("^", subsp), names(chr21.ac))]
      subsp.exome.ac=exome.ac[, grepl(paste0("^", subsp), names(exome.ac))]
      # Calculate chr21 DAF
      subsp.chr21.dac=subsp.chr21.ac[, grepl(".dac$", names(subsp.chr21.ac))]
      subsp.chr21.aac=subsp.chr21.ac[, grepl(".aac$", names(subsp.chr21.ac))]
      subsp.chr21.total.dac=rowSums(subsp.chr21.dac)
      subsp.chr21.total.aac=rowSums(subsp.chr21.aac)
      subsp.chr21.total.daf=subsp.chr21.total.dac/(subsp.chr21.total.dac+subsp.chr21.total.aac)
      # Calculate exome DAF
      subsp.exome.dac=subsp.exome.ac[, grepl(".dac$", names(subsp.exome.ac))]
      subsp.exome.aac=subsp.exome.ac[, grepl(".aac$", names(subsp.exome.ac))]
      subsp.exome.total.dac=rowSums(subsp.exome.dac)
      subsp.exome.total.aac=rowSums(subsp.exome.aac)
      subsp.exome.total.daf=subsp.exome.total.dac/(subsp.exome.total.dac+subsp.exome.total.aac)
      # Remove fixed sites?
      if(fixed.sites==F){
        subsp.exome.total.daf=subsp.exome.total.daf[subsp.exome.total.daf!=0 & subsp.exome.total.daf!=1]
        subsp.chr21.total.daf=subsp.chr21.total.daf[subsp.chr21.total.daf!=0 & subsp.chr21.total.daf!=1]
      }
      # Plot - normal scale
      SFS=ggplot(NULL) + 
        theme_classic() +
        # Density plots (bandwidth set very low to see all the lumps and bumps)
        #geom_density(aes(x=subsp.chr21.total.daf, col='chr21'), alpha=0.3, bw=0.0001) +
        #geom_density(aes(x=subsp.exome.total.daf, col='Exome'), alpha=0.3, bw=0.0001) +
        stat_bin(aes(x=subsp.chr21.total.daf, y=..density.., col='Non-genic\n-chr21'), bins=bins, geom="step", position="identity") +
        stat_bin(aes(x=subsp.exome.total.daf, y=..density.., col='Exome'), bins=bins, geom="step", position="identity") +
        # General
        labs(title=title, subtitle=paste0(subtitle), x="Derived allele frequency (DAF)") +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
              axis.title=element_text(size=15), axis.text = element_text(size=10),
              legend.title=element_text(size=20), legend.text=element_text(size=15), legend.position = c(0.8, 0.9)) +
        scale_discrete_manual(name='', values = c('Exome' = 'turquoise4', 'Non-genic\n-chr21' = 'violetred'), aesthetics = c("colour"))
      print(SFS+labs(y="Density"))
      if(log){print(SFS+labs(y=expression(log[10]~density))+scale_y_log10())}
    }
  }
}

exome.vs.chr21.sfs.genic=function(exome.ac, chr21.ac, fixed.sites=F, title="Exome and chr21 SFS", bins=100, log=TRUE){
  # Get the subspecies prefixes
  chr21_pop_names=names(chr21.ac)
  chr21_pop_names=chr21_pop_names[3:length(chr21_pop_names)]
  chr21_subsps=unique(substr(chr21_pop_names, start = 1 , stop = 1))
  exome_pop_names=names(exome.ac)
  exome_pop_names=exome_pop_names[3:length(exome_pop_names)]
  exome_subsps=unique(substr(exome_pop_names, start = 1 , stop = 1))
  subsps=unique(c(chr21_subsps, exome_subsps))
  if(length(subsps)>1){subsps=c("", subsps)}else{subsps=""}
  for(subsp in subsps){
    if(subsp=='c'){
      subtitle="Central"
    }else if(subsp=='e'){
      subtitle="Eastern"
    }else if(subsp=='n'){
      subtitle="Nigeria-Cameroon"
    }else if(subsp=='w'){
      subtitle="Western"
    }else{subtitle=""}
    if(sum(grepl(paste0("^", subsp), names(chr21.ac)))>1 & sum(grepl(paste0("^", subsp), names(exome.ac)))>1){
      # Select columns corresponding to the subspecies
      subsp.chr21.ac=chr21.ac[, grepl(paste0("^", subsp), names(chr21.ac))]
      subsp.exome.ac=exome.ac[, grepl(paste0("^", subsp), names(exome.ac))]
      # Calculate chr21 DAF
      subsp.chr21.dac=subsp.chr21.ac[, grepl(".dac$", names(subsp.chr21.ac))]
      subsp.chr21.aac=subsp.chr21.ac[, grepl(".aac$", names(subsp.chr21.ac))]
      subsp.chr21.total.dac=rowSums(subsp.chr21.dac)
      subsp.chr21.total.aac=rowSums(subsp.chr21.aac)
      subsp.chr21.total.daf=subsp.chr21.total.dac/(subsp.chr21.total.dac+subsp.chr21.total.aac)
      # Calculate exome DAF
      subsp.exome.dac=subsp.exome.ac[, grepl(".dac$", names(subsp.exome.ac))]
      subsp.exome.aac=subsp.exome.ac[, grepl(".aac$", names(subsp.exome.ac))]
      subsp.exome.total.dac=rowSums(subsp.exome.dac)
      subsp.exome.total.aac=rowSums(subsp.exome.aac)
      subsp.exome.total.daf=subsp.exome.total.dac/(subsp.exome.total.dac+subsp.exome.total.aac)
      # Remove fixed sites?
      if(fixed.sites==F){
        subsp.exome.total.daf=subsp.exome.total.daf[subsp.exome.total.daf!=0 & subsp.exome.total.daf!=1]
        subsp.chr21.total.daf=subsp.chr21.total.daf[subsp.chr21.total.daf!=0 & subsp.chr21.total.daf!=1]
      }
      # Plot - normal scale
      SFS=ggplot(NULL) + 
        theme_classic() +
        # Density plots (bandwidth set very low to see all the lumps and bumps)
        #geom_density(aes(x=subsp.chr21.total.daf, col='chr21'), alpha=0.3, bw=0.0001) +
        #geom_density(aes(x=subsp.exome.total.daf, col='Exome'), alpha=0.3, bw=0.0001) +
        stat_bin(aes(x=subsp.chr21.total.daf, y=..density.., col='chr21'), bins=bins, geom="step", position="identity") +
        stat_bin(aes(x=subsp.exome.total.daf, y=..density.., col='Exome'), bins=bins, geom="step", position="identity") +
        # General
        labs(title=title, subtitle=paste0(subtitle), x="Derived allele frequency (DAF)") +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
              axis.title=element_text(size=15), axis.text = element_text(size=10),
              legend.title=element_text(size=16), legend.text=element_text(size=15), legend.position = c(0.8, 0.85)) +
        scale_discrete_manual(name='Capture', values = c('Exome' = 'turquoise4', 'chr21' = 'violetred'), aesthetics = c("colour"))
      print(SFS+labs(y="Density"))
      if(log){print(SFS+labs(y=expression(log[10]~density))+scale_y_log10())}
    }
  }
}


format_daf=function(df){
  df.dac=as.matrix(df[, grepl(".dac$", names(df))])
  df.aac=as.matrix(df[, grepl(".aac$", names(df))])
  df.daf=as.data.frame(df.dac/(df.dac+df.aac))
  df.daf=cbind(df[,c('chr', 'pos')], df.daf)
  df.daf=melt(setDT(df.daf), id.vars = c("chr","pos"), variable.name = "Population")
  colnames(df.daf)[colnames(df.daf)=='value']="DAF"
  df.daf$chr_pos=paste(df.daf$chr, df.daf$pos, sep="_")
  df.daf=df.daf[order(df.daf$chr_pos),]
  df.daf[df.daf$DAF=="NaN", "DAF"]=NA
  df.daf$Population=gsub(".dac", "", df.daf$Population)
  df.daf$Population=gsub("_", ".", df.daf$Population)
  return(df.daf)
}