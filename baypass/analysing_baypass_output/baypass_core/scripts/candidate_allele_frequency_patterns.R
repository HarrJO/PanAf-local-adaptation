library(dendextend)
library(gplots)
library(cluster)
library(fpc)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(maps)
library(mapdata)
library(mapplots)
library(rgdal)
library(scatterpie)
library(ggrepel)
library(shadowtext)

# Maps

# Downloading and extracting chimp ranges polygons
## Downloading 
### Subspecies ranges were downloaded from https://www.iucnredlist.org/resources/spatial-data-download (downloading the full file for terrestrial mammals). 
### You need to make an account (should stay logged in on google chrome) and request the download. 
### The downloaded file is at `/Users/harrisonostridge/OneDrive - University College London/Projects/PanAf/maps/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp`.
## Extract chimp polygons
### Here I open this file, select only the chimp subspecies and then write this new smaller file. This only needs to be run once. I then read this file in in the next chunk.
#### ranges=readOGR("../../../../maps/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
#### chimp.ranges=ranges[ranges$binomial=="Pan troglodytes",]
#### writeOGR(chimp.ranges, "../../../../maps/chimp_ranges/", "chimp_ranges", driver = "ESRI Shapefile")
#### chimp.ranges=readOGR("../../../../maps/chimp_ranges/chimp_ranges.shp")

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
      geom_scatterpie(data=allele.freq_cand.snp, aes(x=plot.long, y=plot.lat, group=Population, r=0.75), 
                      cols=c("DAF","AAF"), legend_name = "Alleles") +
      geom_text_repel(aes(label=Population), 
                      size = 2, show.legend = FALSE, box.padding=0.25, point.padding=5, nudge_y=1, segment.alpha=0) +
      # General
      labs(title=title, x="", y = "") +
      theme(plot.title = element_text(hjust = 0.5))
    print(plot)
  }
}


# Standradised allele frequencies 
plot_std_af_map=function(allele.freq_cand, chimp.ranges, title_prefix=""){
  # allele.freq_cand requires:
  ## a row per SNP per population
  ## Columns MRK, Population, M_Pstd, Latitude and Longitude
  if(!all(c("MRK", "Population", "M_Pstd", "Latitude", "Longitude") %in% colnames(allele.freq_cand))){
    stop("Dataframe must have columns: MRK, Population, M_Pstd, Latitude and Longitude")}
  # chimp.ranges can be downloaded from https://www.iucnredlist.org/resources/spatial-data-download and chimp information extracted in R with:
  ## ranges=readOGR("../../maps/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
  ## chimp.ranges=ranges[ranges$binomial=="Pan troglodytes",]
  ## writeOGR(chimp.ranges, "../../maps/chimp_ranges/", "chimp_ranges", driver = "ESRI Shapefile")
  ## chimp.ranges=readOGR("../../../../maps/chimp_ranges/chimp_ranges.shp")
  
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
  
  # Point size
  point.size=100/(ymax-ymin)
  # Plot
  world_map=map_data("world")
  
  for(snp in unique(allele.freq_cand$MRK)){
    allele.freq_cand.snp=allele.freq_cand[allele.freq_cand$MRK==snp,]
    ## Title
    if('chr' %in% colnames(allele.freq_cand.snp) & 'pos' %in% colnames(allele.freq_cand.snp)){
      title=paste0(title_prefix, "Standardised Allele Frequencies, chr", unique(allele.freq_cand.snp$chr),":", unique(allele.freq_cand.snp$pos))
    }else{
      title=paste0(title_prefix, "Standardised Allele Frequencies, SNP: ", snp)
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
      # Plot points
      geom_segment(aes(x = Longitude, y = Latitude, xend = plot.long, yend = plot.lat)) +
      geom_point(data=allele.freq_cand.snp, aes(x=plot.long, y=plot.lat, fill=M_Pstd), size=point.size, pch=21, colour='black') +
      scale_fill_gradient2(low = "red", mid = "white", midpoint = 0, high = "green", space = "Lab" ) +
      geom_text_repel(aes(label=Population), 
                      size = 2, segment.color="black", segment.linetype=1, segment.size=0.2, show.legend = FALSE, colour="black") +
      # General
      labs(title=title, x="", y = "") +
      theme(plot.title = element_text(hjust = 0.5))
    print(plot)
  }
}

# Correlation

baypass.2.snp.by.pop=function(allele.freq, stat){
  # Reshape so we have a snp x population matrix 
  snp.by.pop=reshape(allele.freq[,c("MRK","chr", "pos" ,stat, "Population")], idvar = c("MRK","chr", "pos"), timevar = "Population", direction = "wide")
  # Ensure rownames are genopmic coordinates  
  rownames(snp.by.pop)=paste(snp.by.pop$chr, snp.by.pop$pos, sep="_")
  # Remove the stat name from column names
  colnames(snp.by.pop)=gsub(paste0(stat,"."), "", colnames(snp.by.pop))
  # Select only population columns (which start with subspecies prefix)
  pop.cols=grep("^c\\.|^e\\.|^n\\.|^w\\.", colnames(snp.by.pop))
  # Select population columns, convert to matrix and return
  snp.by.pop=as.matrix(snp.by.pop[,pop.cols])
  return(snp.by.pop)
}

make.dendrogram=function(snp.by.pop, method="average", plot=TRUE){
  # Make correlation matrix
  snp.by.pop.cor=cor(snp.by.pop)
  # Make distance matrix
  dist=as.dist(1 - snp.by.pop.cor)
  # Make tree
  tree=hclust(dist, method=method)
  # Make dendrogram
  dend=as.dendrogram(tree)
  # Plot 
  if(plot==TRUE){plot(dend, horiz = TRUE)}
  return(dend)
}

allele.freq.correlation.matrix.4=function(snp.by.pop, title="Correlation", pdf_file=NULL, name.col=NULL, side.col=NULL, cluster_method='none', 
                                          min_k=1, max_k=20, cut_height=0.7){
  if(!(cluster_method %in% c('k_medoids', 'tree_cut','none'))){
    stop("Select either 'k_medoids', 'tree_cut' or none as the clustering method.")
  }
  ## Make dendrogram
  dend=make.dendrogram(snp.by.pop, method="average", plot=FALSE)
  group.col=c()
  if(cluster_method=='k_medoids'){
    # Only works if we have at least 3 columns
    if(ncol(snp.by.pop)>2){
      # pamk - find discrete clusters with K-medoids Clustering
      ## Work out optimal K and find clusters 
      ### pamk documentation: "If 1 is included, a Duda-Hart test is applied and 1 is estimated if this is not significant"
      ### If there are too many categories then pamk takes too long, in this case we make a dissimilarity matrix and run pamk on that instead
      pamk.out=pamk(t(snp.by.pop), krange=min_k:min(max_k, ncol(snp.by.pop)-1))
      ## Assign colours to each cluster
      clust=pamk.out$pamobject[3]$clustering
      colours=colorRampPalette(brewer.pal(name="Paired", n = 12))(max(clust))
      for(pop in clust){group.col=c(group.col, colours[pop])}
    }else{cluster_method='none'}
  } 
  if(cluster_method=='tree_cut'){
    # Select clusters by cutting the dendrogram
    clust=cutree(dend, h=cut_height)
    colours=colorRampPalette(brewer.pal(name="Paired", n = 12))(max(clust))
    for(pop in clust){group.col=c(group.col, colours[pop])}
  } 
  if(!is.null(side.col)){
    group.col=side.col
  }
  if(cluster_method=='none' & is.null(side.col)){
    group.col=rep('white', ncol(snp.by.pop))
  }
  ## Label colours
  if(is.null(name.col) & is.null(name.col)){name.col=rep('black', ncol(snp.by.pop))}
  # Plot heatmap
  ## Required to ensure column and row colours are in the correct order when using revC=T (this makes it much easier to interpret)
  ### from https://stackoverflow.com/questions/34136096/r-gplots-heatmap-with-side-colours
  ### First make a heatmap without reversing column order
  pdf(file = NULL) ### This ensures the plot is not displayed 
  dummy.heatmap=heatmap.2(cor(snp.by.pop), Rowv = ladderize(dend), Colv = ladderize(dend), revC = FALSE,dendrogram = "both", RowSideColors=group.col)
  dev.off()
  #### Get row order
  ordinary_order = dummy.heatmap$rowInd
  #### Reverse row order
  reversal = cbind(ordinary_order, rev(ordinary_order))
  ##### Column/row colours
  group.col_rev = group.col[reversal[,2]]
  group.col_rev = group.col_rev[order(reversal[,1])]
  ##### Column/row labels
  labs_rev = colnames(snp.by.pop)[reversal[,2]]
  labs_rev = labs_rev[order(reversal[,1])]
  ##### Column/row label colours
  name.col_rev = name.col[reversal[,2]]
  name.col_rev = name.col_rev[order(reversal[,1])]
  ## Plot actual heatmap
  #col = colorRampPalette(c("blue", "yellow", "red"))(500)
  col = plasma(100)
  text.size=min(1.5, 1/(ncol(snp.by.pop)*0.07))
  if(!is.null(pdf_file)){pdf(file = pdf_file)}
  cex=min(20/nrow(cor(snp.by.pop)), 1.5)
  par(oma=c(0,0,1.5,0))
  heatmap.2(cor(snp.by.pop),
            Rowv = ladderize(dend),  
            Colv = ladderize(dend), 
            revC = TRUE, # This means the RowSideColors needs to be reversed to stay in sync
            dendrogram = "both", 
            col = col,
            trace = "none", density.info = "density",
            margins = c(15, 15),
            key.title=NA, key.ylab=NA, key.xlab=NA, keysize=1,
            main=title,
            ColSideColors=group.col, RowSideColors=group.col_rev,
            cexRow=cex, cexCol=cex,
            colRow=name.col_rev, colCol=name.col, 
            labRow=labs_rev, labCol=colnames(snp.by.pop)
  )
  if(!is.null(pdf_file)){dev.off()}
  if(cluster_method!='none'){
    return(list(dendrogram=dend, clusters=clust))
  }
}

# Correlate standardised allele frequencies with geography
## Here I use PC axes to describe geography and test for correlation between M_Pstd and these PCs
cor_geography.vs.AF=function(allele.freq, title_prefix="", chimp.ranges){
  # Select geographical coordinates
  coord=unique(allele.freq[,c('Population', 'Longitude', 'Latitude')])
  rownames(coord)=coord$Population
  coord=as.matrix(coord[,c('Longitude', 'Latitude')])
  # Run PCA
  ## We do this so we can use PC1 as a proxy of geography so we can test for a correlation between allele frequencies and geography
  coord_pca=prcomp(t(coord), scale = FALSE)
  ## Get variation per PC
  v = coord_pca$sdev^2
  pv = 100*v/sum(v)
  ## Add PC coordinates to data frame
  coord_pca_r=as.data.frame(coord_pca$rotation)
  coord_pca_r$Population=rownames(coord_pca_r)
  # Add PC coordinates
  tmp=merge(allele.freq, coord_pca_r, all.x=TRUE, by='Population')
  # For each PC
  for(PC in 1:2){
    # Test for correlation
    cor_test_res=cor.test(tmp[['M_Pstd']], tmp[[paste0("PC", PC)]])
    # Plot correlation 
    plot=ggplot(tmp, aes_string(x=paste0("PC", PC), y='M_Pstd')) +
      theme_bw() +
      geom_point(aes(col=Population)) +
      geom_smooth(formula=y~x, method='lm', col="black") +
      labs(title=paste0(title_prefix, "Correlation coefficent: ", round(cor_test_res$estimate,3), ",   p-value: ", round(cor_test_res$p.value,3)),
           x=paste0("PC", PC, " of geographical coordinates (",pv[PC],"% variation)")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(plot)
  }
}

# Plot PCA
af_PCA=function(snp.by.pop, pop.info, PCs=c(1,2), colour="Habitat", title=""){
  # Closely followed https://www.youtube.com/watch?v=0Jp4gsfOLMs
  # Performing PCA
  ## prcomp expects tables to have samples as rows and SNPs as columns so I transpose the data frame 
  ## prcomp manual says that while scaling variables to have unit variance (scale = TRUE) is advisable, it cannot be used if there are zero values (which I have). Scaling the variables results in unusual PCAs.
  p = prcomp(t(snp.by.pop), scale = FALSE)
  v = p$sdev^2
  pv = 100*v/sum(v)
  # Scree plot (Bar plot)
  pc.var=data_frame("PC"=1:length(pv), "variance"=pv)
  pc.var$rounded_variance=round(pc.var$variance, 2)
  ## Ensure bars will be in the correct order
  pc.var$PC=factor(pc.var$PC, levels = pc.var$PC)
  bar=ggplot(pc.var, aes(x=PC, y=variance)) + 
    ## General
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size=7),
      axis.text.x = element_text(size=7)) +
    labs(title=paste0("% Variance Explained by Each PC: ", title), x="PC", y="Percentage Variance") +
    ## Plot
    geom_bar(stat = "identity") +
    # Text
    geom_text(aes(label=rounded_variance), vjust=-1, color="black", size=2)
  print(bar)
  # Make PCA data frame for plotting
  pca.df=data.frame(p$x)
  pca.df$Population=rownames(pca.df)
  pop.info=unique(pop.info[,c('Population', colour)])
  pca.df=merge(pca.df, pop.info)
  # Plot
  plot=ggplot(pca.df, aes_string(x=paste0("PC", PCs[1]), y=paste0("PC", PCs[2]), fill=colour, colour=colour, label = "Population")) +
    # General
    theme_minimal()+ 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.key.size = unit(0.4, "cm"), 
          plot.title = element_text(hjust = 0.5),
          legend.position="none") +
    # Plot
    geom_point(alpha=0.1) +
    geom_text_repel(size = 3, min.segment.length = unit(0.1, "lines"), show.legend = FALSE) +
    # Labels and legends 
    labs(title=paste0("PCA Plot: ", title))
  if(colour %in% c("Habitat", "Subspecies")){
    plot=plot+scale_fill_manual(values = c('Central'='green3', 'Eastern'='darkorange', 'Nigeria-Cameroon'='red', 'Western'='blue',
                                           'forest'='darkgreen', 'mosaic'='darkgoldenrod4', 'savanna'='khaki3')) +
      scale_color_manual(values = c('Central'='green3', 'Eastern'='darkorange', 'Nigeria-Cameroon'='red', 'Western'='blue',
                                    'forest'='darkgreen', 'mosaic'='darkgoldenrod4', 'savanna'='khaki3'))
  }
  print(plot)
}

af_PCA_each_tail=function(allele.freq, xtx, tails=c(1, 0.005, 0.001, 0.0005), stat, PCs=c(1,2)){
  for(tail in tails){
    if(tail==1){
      allele.freq.tail=allele.freq
    }else{
      allele.freq.tail=allele.freq[allele.freq$MRK %in% xtx[xtx$empirical_p<tail,]$MRK,]
    }
    snp.by.pop.tail=baypass.2.snp.by.pop(allele.freq.tail, stat)
    af_PCA(snp.by.pop.tail, pop.info=allele.freq,PCs=PCs, colour="Subspecies", title="Test")
  }
}


# Plot PCA biplot
PCbiplot=function(PC, PCs=c(1,2), title="", point.cols='black', var_bar=FALSE, min_label_size=1) {
  
  # From https://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
  # PC being a prcomp object
  
  # Scree plot (Bar plot)
  ## Get variance explained by PCs
  v = PC$sdev^2
  pv = 100*v/sum(v)
  pc.var=data_frame("PC"=1:length(pv), "variance"=pv)
  pc.var$rounded_variance=round(pc.var$variance, 2)
  ## Ensure bars will be in the correct order
  pc.var$PC=factor(pc.var$PC, levels = pc.var$PC)
  if(var_bar){
    bar=ggplot(pc.var, aes(x=PC, y=variance)) + 
      ## General
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=7),
        axis.text.x = element_text(size=7)) +
      labs(title=paste0("% Variance Explained by Each PC: ", title), x="PC", y="Percentage Variance") +
      ## Plot
      geom_bar(stat = "identity") +
      # Text
      geom_text(aes(label=rounded_variance), vjust=-1, color="black", size=2)
    print(bar)
  }
  # PC plot
  x=paste0("PC", PCs[1])
  y=paste0("PC", PCs[2])
  data=data.frame(obsnames=row.names(PC$x), PC$x)
  ## Plot
  plot=ggplot(data, aes_string(x=x, y=y)) + 
    theme_minimal()+ 
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") + 
    geom_hline(yintercept=0, size=0.2) + geom_vline(xintercept=0, size=0.2) +
    labs(title=title, x=paste0(x, " (", round(pv[PCs[1]], 2) ,"%)"), y=paste0(y, " (", round(pv[PCs[2]], 2) ,"%)"))
  ### Arrows
  datapc=data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult=min((max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
           (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x]))))
  #### Control the length of arrows (I think the length is arbitrary - all that matters is relative length)
  length_multiplyer=1.5
  datapc=transform(datapc,v1 = length_multiplyer* mult * (get(x)),v2 = length_multiplyer* mult * (get(y)))
  #### Add arrows
  arrow_transpareny=min(1/(0.01*nrow(datapc)), 0.3)
  plot=plot + coord_equal() + 
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=arrow_transpareny, color="red")
  ### Arrow and point labels
  arrow_label_size=min(50/nrow(datapc), 3)
  label_size=min(90/nrow(data), 3)
  if(label_size<min_label_size){
    #plot=plot+geom_point(alpha=0.75, shape=16, size=0.5, col=point.cols)
    plot=plot+geom_point(alpha=1, shape=16, size=1, col=point.cols)
  }else{
    #### geom_shadowtext() gives a black outline around the text so you can see it even if in a light colour such as yellow
    plot=plot+geom_shadowtext(alpha=1, size=label_size, aes(label=obsnames), col=point.cols)
  }
  if(arrow_label_size>min_label_size){
    plot=plot+geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = arrow_label_size, vjust=1, color="red")
  }
  # Print plot
  print(plot)
}

