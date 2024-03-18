# simplified molecule distance script 
# calculate the distance between two labels in the specified input files. 
# optimized to handle Guided Assembly & RVA outputs (alignmolvref/merge/) 
# output a plot of the distances

#' @param rcmappath : anchor/reference cmap
#' @param xmap: xmap with alignment information
#' @param qcmappath: query/molecule cmap
#' @param startlabelid: ID of left reference label of interest
#' @param endlabelid: ID of right reference label of interest
#' @param out_handle: prefix for the output files
#' @return A table of distances between labels IDs (out_handle.csv)


# 
library(ggplot2)

# command line args
cmd_args <- commandArgs(trailingOnly=TRUE)
rcmappath <- cmd_args[1]
xmap <- cmd_args[2]
qcmappath <- cmd_args[3]
startlabelid <- cmd_args[4]
endlabelid <- cmd_args[5]
out_handle <- cmd_args[6]

##
# helper functions
##
correspondingmolcoor<-function(rownum, reflabelnum, withalignments_1){
  permolalignment<-withalignments_1[rownum,]
  eachlabel<-strsplit(permolalignment$Alignment, split="[)(]")
  wantedqrylabels<-eachlabel
  
  getperwantedqrylabel<-function(num, reflabelnum){
    perwantedqrylabels<-wantedqrylabels[[num]]
    the_align<-perwantedqrylabels[grep(reflabelnum, perwantedqrylabels)]
    the_align<-unlist(strsplit(the_align, split=". "), recursive=F)
    peralign<-strsplit(the_align, split=",")
    findcorrespond<-data.frame(matrix(unlist(peralign), length(peralign), byrow =T))
    good<-which(findcorrespond$X1 == reflabelnum)
    if(length(good)>0){
      the_align<-the_align[good]
      qrylabel<-do.call(rbind, strsplit(the_align, split=","))[,2]
      
      return(qrylabel)
    } else {
      return("NA")
    }
  }
  qrylabels<-lapply(seq(1, length(wantedqrylabels), 1), getperwantedqrylabel, reflabelnum)
  return(qrylabels)
}

getqrycoor<-function(rcmap, xmap, labelstart){

  r<-rcmap
  x<-xmap
  
  refstart<-r[r$SiteID == labelstart,]
  refstartpos<-data.frame(as.numeric(as.character(refstart$CMapId)), refstart$Position)
  colnames(refstartpos)<-c("CMapId", "Position")
  
  refstartmolcount<-refstart$Occurrence
  withalignments<-x[grepl(labelstart, x$Alignment),]
  print("var initialized")
  
  # print(str(refstart))
  # print(str(refstartpos))

  num_row<-as.integer(nrow(withalignments))  
  qrycoor<-lapply(seq(1, num_row, 1), correspondingmolcoor, labelstart, withalignments_1=withalignments)
  
  result<-data.frame(withalignments$RefContigID ,withalignments$QryContigID, do.call(rbind, qrycoor))
  colnames(result)<-c("RefContigID", "QryContigID", "Qrylabel")
  print("result data frame initiated")
  print(str(result))
  
  result$refPos<-refstartpos$Position[match(result$RefContigID, refstartpos$CMapId)]
  # result<-merge(result, refstartpos, by.x )
  print("refPos update")
  print(str(result))
  
  result$Qrylabel<-unlist(result$Qrylabel)
  colnames(result)<-c("RefContigID", "QryContigID", "Qrylabel", "refpos")  
  result<-result[c("QryContigID", "Qrylabel", "refpos")]
  result$Qrylabel<-unlist(result$Qrylabel)
  print("result subset")
  print(str(result))
  
  result<-result[!result$Qrylabel=="NA",]
  print(str(result))
  
  return(result)
}

getcoor<-function(n, starts, qcmap_1){
  perqry<-starts[n,]
  findcoor<-qcmap_1[qcmap_1$CMapId == perqry$QryContigID,]
  findcoor<-findcoor[findcoor$SiteID == perqry$Qrylabel,]$Position
  return(findcoor)
}

getcolnames<-function(rarepath){
  rare1<-read.table(rarepath, comment.char="#", sep="\t", stringsAsFactors = F)
  nhl<-system(paste0("grep '^#' ", rarepath, "| wc -l"), intern =TRUE)
  header<-system(paste0("less ", rarepath, "| grep '^#h'"), intern=TRUE)
  if(length(header)==0){
    header<-system(paste0("less ", rarepath, "| grep '^# SVIndex'"), intern=TRUE)
  } 
  if(length(header)==0){
    header<-system(paste0("less ", rarepath, "| grep '^#clusterId'"), intern=TRUE)
  } 
  if(length(header)==0){
    header<-system(paste0("less ", rarepath, "| grep '^#'"), intern=TRUE)
  } 
  colnames(rare1)<-strsplit(header, split="\t")[[1]]
  return(rare1)
}

getheader<-function(rarepath){
  rare1<-read.table(rarepath, comment.char="#", sep="\t", stringsAsFactors = F)
  nhl<-system(paste0("grep '^#' ", rarepath, " | wc -l"), intern =TRUE)
  header<-readLines(rarepath, n=nhl)
  return(header)
}

waterfall_plotting_distance <- function(distance_data, filename_prefix){
  distance_mean <- mean(distance_data$distance)
  distance_data$index <- row.names(distance_data)
  waterfall_plot<- ggplot(distance_data, aes(x=reorder(index,distance), y=distance)) +
    geom_bar(stat="identity", 
             fill="#0072B2",
             width=0.7,
             position = position_dodge(width=0.4))+
    scale_y_discrete(expand = c(0,0)) +
    labs(title = "Distance Between Specified Labels", 
         x = "Molecules", 
         y = "Distance (basepairs)") +
    theme_classic()+
    theme(axis.text.x = element_blank()) +
    theme(legend.position="none") +
    geom_hline(yintercept=distance_mean, linetype="dashed", color = "#D55E00")+
    geom_text(aes(20,distance_mean,label = "mean", vjust = -1), 
             color = "#D55E00", check_overlap=T )
    ggsave(plot = waterfall_plot, file = paste0(filename_prefix,'_barplot.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300)
  }

waterfall_plotting_delta <- function(distance_data, filename_prefix){
  distance_data$index <- row.names(distance_data)
  distance_mean <- mean(distance_data$delta)
  waterfall_plot<- ggplot(distance_data, aes(x=reorder(index,delta), y=delta)) +
    geom_bar(stat="identity", 
             fill="#0072B2",
             width=0.7,
             position = position_dodge(width=0.4))+
    # scale_y_discrete(expand = c(0,0)) +
    labs(title = "Molecule Label Differences From hg38 Reference", 
         x = "Molecules", 
         y = "Δ Distance (Molecule - Reference) (basepairs)") +
    theme_classic()+
    theme(axis.text.x = element_blank()) +
    theme(legend.position="none") +
    geom_hline(yintercept=distance_mean, linetype="dashed", color = "#D55E00")+
    geom_text(aes(20,distance_mean,label = "mean", vjust = -1), color = "#D55E00", check_overlap=T )
    ggsave(plot = waterfall_plot, file = paste0(filename_prefix,'_delta_barplot.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300)
  }

summary_stats <- function(distance_data, filename_prefix){
  mean_distance<-mean(distance_data$distance)
  sd_distances<- sd(distance_data$distance)
  cv_distances <- sd_distances/mean_distance
  summary_table <- cbind(mean_distance, sd_distances,cv_distances)
  write.csv(file=paste0(filename_prefix, '_summary_stats.csv'), x=summary_table)
  print(summary_table)
}

histogram_plot <- function(distance_data, filename_prefix){
  num_bins <- sqrt(nrow(distance_data))
  hist_dens_plot <- ggplot(distance_data, aes(x = distance)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "bisque1", bins=num_bins) +
    geom_density()+
    labs(title = "Histogram of Distances with Density Estimate", 
         x = "Distance (basepairs)", 
         y = "Density") +
    theme_classic()
    ggsave(plot = hist_dens_plot, file = paste0(filename_prefix,'_histogram.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300) 
}

violin_plot <- function(distance_data, filename_prefix){
  distance_data$arbitrary_group <- 'xyz'
  violin_with_boxplot <- ggplot(distance_data, aes(y=distance, x=arbitrary_group)) + 
    geom_violin(trim=FALSE, fill="darkslategray4")+
    labs(title="Violin Plot of Distances (with Boxplot)", x = "Sample", y = "Distances (basepairs)")+
    geom_boxplot(width=0.1)+
    theme_classic()+ 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
    ggsave(plot = violin_with_boxplot, file = paste0(filename_prefix,'_violin_plot.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300)
}

clustering <- function(completeframe, 
                        startlabelid, 
                        endlabelid, 
                        out_handle,nk, 
                        additional_noise_removal=FALSE, 
                        binw=400){
    completeframeoriginal<-round(completeframe$distance,0)
    
    mod <- Mclust(completeframe$distance, modelNames="V")
    if(!is.null(nk)){
      mod<-Mclust(completeframe$distance, modelNames="V", G=nk)
    }
    j<-summary(mod)
    summary(mod$BIC)
    mclust1Dplot(data = completeframe$distance, what = "density", parameters=mod$parameters, z=mod$z)
    abline(v=mod$parameters$mean, lty=3)
    mclust1Dplot(data = completeframe$distance, what = "classification", parameters=mod$parameters, z=mod$z)
    mclust1Dplot(data = completeframe$distance, what = "uncertainty", parameters=mod$parameters, z=mod$z)  ###uncertainty is how uncertain it is to be asscoiated with the group classified. 
    plot.mclustBIC(mod$BIC)
    k<-mod$G
    #print(nrow(completeframe))
    #print(paste(completeframe$refdistance[1], "ref distance"))
    print(mod$parameters$mean)
    print(mod$parameters$pro)
    print(mod$sd)
    print(mod$parameters$variance$sigmasq)
    hist(mod$uncertainty)
    print("done round 1")
    colors<-c("red", "green", "blue")
    
    
    ###Clear data in outlier cluster or reduce number of clusters
    dis<-data.frame(table(mod$classification))
    
    ###Plotting function
    Plotting_hist_and_Gau<-function(mod, startlabelid, endlabelid){
      newcolor<-array()
     # print(colors)
      adjusteddnorm<-function(data, mm, sd, prop){
        re<-prop*dnorm(data, mean=mm, sd=sd)*nrow(completeframe)*binw
        return(re)
      }
      
      p1<-ggplot(data.frame(completeframeoriginal), aes(x=completeframeoriginal))+geom_histogram(aes(x=completeframeoriginal, y=..count..), fill=NA, color="black", binwidth=binw)
      
      #print(length(mod$parameters$mean))
      
      ###Gather Component info for legend
      legwrite=array()
      for (k in 1:length(mod$parameters$mean)){
        #print(k)
        newleg<-paste0("C-", nrow(completeframe), "\nD-", completeframe$refdistance[1], "\nM-", format(mod$parameters$mean[k],digits=3), "\nSD-", format(mod$sd[k], digits=3), "\nSD/M-", format((mod$sd[k]/mod$parameters$mean[k]*100), digit=3), "%\nP-", format(mod$parameters$pro[k], digits=2))
        #print(newleg)
        legwrite<-append(legwrite, newleg)
      }
      legwrite=legwrite[!is.na(legwrite)]
      legwrite=rev(legwrite)
      
      for (i in 1:length(mod$parameters$mean)){
        print(paste("loop", i))
        print(mod$G)
        
        #print(colors[i])
        #legwrite=c("G1", "G21")
        #print(legwrite)
        newcolor<-append(newcolor, colors[i])
        
        p1<-p1+
          # stat_function(fun=adjusteddnorm, n=1000, args=list(mm=mod$parameters$mean[i], sd=mod$sd[i], prop=mod$parameters$pro[i]), aes(fill=!!colors[i], color=!!colors[i]), geom="area", alpha = .2)+
          ggtitle(paste0("clustering molecule distances between", "\n labels ", startlabelid, "-", endlabelid))+
          xlab("Distance (bp) between Labels of Interest")+scale_color_discrete(guide=FALSE)+ylab("Number of Molecules")+
          geom_vline(xintercept=27686, linetype="dashed", color = "#D55E00")+
          xlim(25000,30000)+
          theme(plot.title = element_text(hjust=0.5)) + theme(legend.position="none")
      } 
      
      print("do plots")
      #print(p1)
      
      for (j in 2:length(ggplot_build(p1)$data)){
        colori<-ggplot_build(p1)$data[[j]]$colour
        #print(j)
        #print(unique(colori))
        newcolor<-unique(append(newcolor, colori))
      }
      newcolor<-newcolor[!is.na(newcolor)]
      #print(p1+scale_fill_manual("Components\nInfo", values=newcolor, labels=legwrite)+scale_color_manual(values=newcolor, labels=legwrite, guide=FALSE))
      ggsave(paste0(out_handle,'_auto_clustering.png'), p1, device='png', dpi=300)
	#return(p1)
    }
    
    while (TRUE %in% (dis$Freq < round(0.15*sum(dis$Freq),0))) {
      ###If one of the cluster is with <15% of mol
      print(paste(k, "clusters were found and needs to be reduced"))
      Gtoberemoved<-which((dis$Freq < 0.15*sum(dis$Freq)) ==TRUE)
      kremove<-dis[dis$Freq < 0.15*sum(dis$Freq),]$Var1
      k<-k-length(kremove)
      if (TRUE %in% ((which(mod$parameters$variance$sigmasq==max(mod$parameters$variance$sigmasq))) == Gtoberemoved)){
        ###If the one cluster having <15% of mol is the cluster with the largest variance
        ###Then remove data of that cluster
        filtereddata<-completeframe[!mod$classification %in% Gtoberemoved,]
      } else {
        filtereddata<-completeframe
      }
      
      #filtereddata<-completeframe
      print(paste(nrow(completeframe)-nrow(filtereddata), "data is removed"))
      mod <- Mclust(filtereddata$distance, modelNames="V", G=k)
      j<-summary(mod)
      summary(mod$BIC)
      mclust1Dplot(data = filtereddata$distance, what = "density", parameters=mod$parameters, z=mod$z)
      abline(v=mod$parameters$mean, lty=3)
      mclust1Dplot(data = filtereddata$distance, what = "classification", parameters=mod$parameters, z=mod$z)
      #mclust1Dplot(data = completeframe$distance, what = "uncertainty", parameters=mod$parameters, z=mod$z)  ###uncertainty is how uncertain it is to be asscoiated with the group classified. 
      plot.mclustBIC(mod$BIC)
      mod$sd<-sqrt(mod$parameters$variance$sigmasq)
      #print(nrow(filtereddata))
      #print(paste(completeframe$refdistance[1], "ref distance"))
      #print(mod$parameters$mean)
      #print(mod$parameters$pro)
      #print(mod$sd)
      #print(mod$parameters$variance$sigmasq)
      #hist(mod$uncertainty)
      dis<-data.frame(table(mod$classification))
      completeframe<-filtereddata
      Plotting_hist_and_Gau(mod, startlabelid, endlabelid)
    } 
    ###Plotting begins here
    
    
    if ((TRUE %in% (diff(sort(mod$parameters$mean, decreasing=FALSE)) <max(sqrt(mod$parameters$variance$sigmasq)))) & (TRUE %in% (dis$Freq < 0.10*sum(dis$Freq)))){
      ###If the difference between means is smaller than the max of the std dev & that one of the clusters is less than 10% of molecules. 
      
      print(paste(k, "clusters were found and needs to be reduced"))
      k<-k-1
      filtereddata<-completeframe
      mod <- Mclust(filtereddata$distance, modelNames="V", G=k)
      j<-summary(mod)
      summary(mod$BIC)
      mclust1Dplot(data = filtereddata$distance, what = "density", parameters=mod$parameters, z=mod$z)
      abline(v=mod$parameters$mean, lty=3)
      mclust1Dplot(data = filtereddata$distance, what = "classification", parameters=mod$parameters, z=mod$z)
      
      plot.mclustBIC(mod$BIC)
      mod$sd<-sqrt(mod$parameters$variance$sigmasq)
      mod$nrow<-nrow(mod$parameters)
      print(nrow(filtereddata))
      print(paste(completeframe$refdistance[1], "ref distance"))
      print(mod$parameters$mean)
      print(mod$parameters$pro)
      print(mod$sd)
      print(mod$parameters$variance$sigmasq)
      hist(mod$uncertainty)
      dis<-data.frame(table(mod$classification))
      completeframe<-filtereddata
      Plotting_hist_and_Gau(mod, startlabelid, endlabelid)
    } else {
      print("I am inside 1")
      mod$sd<-sqrt(mod$parameters$variance$sigmasq)
      ###Plotting begins here
      Plotting_hist_and_Gau(mod, startlabelid, endlabelid)
    }
 #   summary_table <- summary(mod)
    result<-c(startlabelid, completeframe$refdistance[1], length(mod$data), k, mod$parameters$mean, mod$sd, mod$parameters$pro)
    result_summary<-data.frame(c(mod$parameters$mean, mod$sd))
    #print('mean')
    #print(mod$parameters$mean)
    #print('sd')
    #print(mod$sd)
    simple_result<-data.frame(cbind(mod$parameters$mean, mod$sd))
    colnames(simple_result)<-c('mean','sd')
    simple_result$cv <- simple_result$sd/simple_result$mean
    simple_result$total <- mod$n
    simple_result$n <- mod$parameters$pro*mod$n
    print(simple_result)
    write.csv(simple_result, paste0(out_handle, "auto_clustering_cluster_simple_result_.csv"))
    # completeframe$ClusterNum <- mod$classification
    # cluster_summary <-completeframe[,c("ClusterNum", "distance")]
    # write.csv(mod,
    #                 file=paste0(out_handle, "mod_cluster_summary_noOutlier_delta.csv"),
    #           row.names=FALSE,
    #           quote=FALSE)
    
    return(result)
}

### "main" function
molecule_distance_main<-function(rcmappath,
                                 xmap,
                                 qcmappath,
                                 startlabelid, 
                                 endlabelid, 
                                 out_handle
                                 #additional_noise_removal=FALSE # is this needed bc I removed it; was in Joyce's script but we've always had it "off" (e.g. set to FALSE)
                                 )
{
  print("INPUT FILES:")
  print(rcmappath)
  print(xmap)
  print(qcmappath)
  qcmap<-getcolnames(qcmappath)
  qcmap$CMapId<-as.character(qcmap$`#h CMapId`)
  print(rcmappath)
  r<-getcolnames(rcmappath)
  r$CMapId<-as.character(r$`#h CMapId`)
  x<-read.table(xmap, stringsAsFactors=F)
  h<-getheader(xmap)
  h2<-h[length(h)-1]
  h<-strsplit(h2, split="\t")[[1]]
  colnames(x)<-h
  
  starts<-getqrycoor(r, x, startlabelid)

  print("starts done")
  print(str(starts))

  ends<-getqrycoor(r, x, endlabelid)

  print("ends done")
  print(str(ends))
  
  ### Downsample to the most 1000 molecules
  if(dim(starts)[1]>1000){
    starts<-starts[sample(nrow(starts),1000),]
    ends<-ends[sample(nrow(ends), 1000),]
  }
  
  ### get the righmost label if the mols that have multiple alignments to one ref label
  starts$Qrylabel<-unlist(lapply(starts$Qrylabel, max))
  ### get the leftmost label if the mols have multiple alignments to one ref label
  ends$Qrylabel<-unlist(lapply(ends$Qrylabel, min))
  
  is.error <- function(x) inherits(x, "try-error")
  
  ### Get Qry position
  if(nrow(starts)<10 | nrow(ends)<10){   ###Too few molecules at the start of end (Added on 04/20/2021)
    result<-c(r$CMapId[1], startlabelid, starts$refpos[1], min(nrow(starts), nrow(ends)), 0, 0, 0, 0)
    print(paste("Only <10 molecules are the label start or label end"))
    return (result)
  }
  
  starts$QryPos<-lapply(seq(1, nrow(starts), 1), getcoor, starts, qcmap_1=qcmap)
  ends$QryPos<-lapply(seq(1, nrow(ends),1), getcoor, ends, qcmap_1=qcmap)
  completeframe <- merge(starts,ends,by="QryContigID")
  print(paste0(nrow(starts)," molecules at start, ", nrow(ends), " molecules at end, ", nrow (completeframe), " spanned both labels"))
  
  ### x is left breakpoint, y is right breakpoint
  extramol<-which(completeframe$QryPos.x == "numeric(0)"|completeframe$QryPos.y == "numeric(0)")  #Remove any mol that is not in q.cmap (caused by potential RefAligner bug)
  if(length(extramol)>0){
    completeframe<-completeframe[-extramol,]
  }
  completeframe$QryPos.x<-unlist(completeframe$QryPos.x)
  completeframe$QryPos.y<-unlist(completeframe$QryPos.y)
  
  if(nrow(completeframe)<10){   ###Too few molecules to be meaningful
    k<-0
    kmax<-0
    result<-c(r$CMapId[1], startlabelid, completeframe$refdistance[1], nrow(completeframe), k, rep(0, kmax), rep(0, kmax), rep(0, kmax))
    print(paste("Only", nrow(completeframe), " molecules span the labels"))
    return (result)
  } else if(nrow(completeframe)>10){
    
    completeframe$distance<-abs(as.numeric(as.character(completeframe$QryPos.y))-as.numeric(as.character(completeframe$QryPos.x)))
    completeframe$refstartlabelid<-startlabelid
    completeframe$refendlabelid<-endlabelid
    completeframe$refdistance<-abs(as.numeric(as.character(completeframe$refpos.y))-as.numeric(as.character(completeframe$refpos.x)))
    
    completeframe$distance<-round(completeframe$distance,0)

    # get delta of molecule distances from ref distances
    completeframe$delta<-completeframe$distance-completeframe$refdistance
    completeframeoriginal<-round(completeframe$distance,0)
    simple_complete_output <- completeframe[c("QryContigID","distance")]
    write.csv(file = paste0(out_handle, "_distances.csv"), x = simple_complete_output)
    write.csv(file = paste0(out_handle, "_complete_data.csv"), x = completeframe)
    print("QC pass")
  }
  
  waterfall_plotting(simple_output, out_handle)
  histogram_plot(simple_output, out_handle)
  violin_plot(simple_output, out_handle)
  summary_stats(simple_output, out_handle)

}

# function call
molecule_distance_main(rcmappath, xmap, qcmappath, startlabelid, endlabelid, out_handle)

