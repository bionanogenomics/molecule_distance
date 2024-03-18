# Calculate and plot delta of molecule cmap to reference
# Finds delta of molecules to reference sequence and estimates # of repeats
# Outputs multiple plots on molecule distances

#' @param data: table of molecules containing repeats (complete_data.csv)
#' @param startlabelid: ID of left reference label of interest
#' @param endlabelid: ID of right reference label of interest
#' @param ref_length: bp length of reference genome (hg38_DLE1_0kb_0labels.cmap) 
#' @param repeat_length: bp length of repeats (e.g GAA = 3, TCGA = 4)
#' @param out_handle: prefix for the output files
#' @return A plot of distances between labels IDs (delta_out_handle.png)
#' @return A plot of distances between labels IDs (delta_out_handle.csv)

library(mclust)
library(prabclus)
library(ggplot2)

# command line args
cmd_args <- commandArgs(trailingOnly=TRUE)
input_filepath <- cmd_args[1]
startlabelid <- cmd_args[2]
endlabelid <- cmd_args[3]
out_handle <- cmd_args[4]

# main
molecule_distance_delta_ref_main<-function(input_filepath, startlabelid, endlabelid, out_handle){
  # open .csv
  print(out_handle)
  data<-read.csv(input_filepath)
  
  # make new column 'delta' = distance - refdistance
  data$delta<-data$distance-data$refdistance
  
  waterfall_plotting_delta(data, paste0(out_handle,"_complete"))

  # keep deltas within 3 sd
  print(delta_mean<-mean(data$delta))
  print(delta_sd<-sd(data$delta))
  data<-subset(data, abs(data$delta) < (delta_mean+(3*delta_sd)))

  # data<-remove_sd_outlier(data, cols='delta', n_sigmas=3)

  # make new column 'num_repeats' for only +ve distances
  # data$num_repeats<-ifelse(data$delta>0, data$delta/as.integer(repeat_length), 0)
  # data$num_repeats<-data$delta/as.integer(repeat_length)

  # validating data
  print(head(data))
  print(data$num_repeats)
  print(mean(data$distance))

  waterfall_plotting_distance(data, paste0(out_handle,"_trimmed"))
  waterfall_plotting_delta(data, paste0(out_handle,"_trimmed"))
  # waterfall_plotting_numrepeats(data, paste(out_handle,"no_outliers"))

  # Output data
  write.csv(file = paste0(out_handle, "_trimmed_data.csv"), x = data)
  clustering(data,
              startlabelid, 
              endlabelid, 
              out_handle,
              nk = NULL, 
              additional_noise_removal=FALSE, 
              binw=400)

}


# plot histogram of deltas, with line at ref distance
#   x-axis: molecule (shortest -> longest)
#   y-axis: distance from ref
#   bar label: # of repeat (delta/repeat_unit)
waterfall_plotting_distance <- function(distance_data, filename_prefix){
  distance_data$index <- row.names(distance_data)
  distance_mean <- mean(distance_data$distance)
  waterfall_plot<- ggplot(distance_data, aes(x=reorder(index,distance), y=distance)) +
    geom_bar(stat="identity", 
             fill="#0072B2",
             width=0.7,
             position = position_dodge(width=0.4))+
    scale_y_continuous(breaks = round(seq(0, max(distance_data$distance), by = 1000),0))+
    labs(title = "Distance Between Specified Labels", 
         x = "Molecules", 
         y = "Distance (basepairs)") +
    theme_classic()+
    theme(axis.text.x = element_blank()) +
    theme(legend.position="none") # +
    # geom_hline(yintercept=distance_mean, linetype="dashed", color = "#D55E00")+
    # geom_text(aes(10,distance_mean,label = "mean", vjust = -1), 
            #  color = "#D55E00", check_overlap=T )
    ggsave(plot = waterfall_plot, file = paste0(filename_prefix,'_distance_barplot.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300)
  }

waterfall_plotting_delta <- function(distance_data, filename_prefix){
  distance_data$index <- row.names(distance_data)
  distance_mean <- mean(distance_data$delta)
  waterfall_plot<- ggplot(distance_data, aes(x=reorder(index,delta), y=delta)) +
    geom_bar(stat="identity", 
             fill="#0072B2",
             width=0.7,
             position = position_dodge(width=0.4))+
    scale_y_continuous(breaks = round(seq(min(distance_data$delta), max(distance_data$delta), by = 1000),0))+
    labs(title = "Molecule Label Differences From hg38 Reference", 
         x = "Molecules", 
         y = "Δ Distance (Molecule - Reference) (basepairs)") +
    theme_classic()+
    theme(axis.text.x = element_blank()) +
    theme(legend.position="none") # +
    # geom_hline(yintercept=distance_mean, linetype="dashed", color = "#D55E00")+
    # geom_text(aes(10,distance_mean,label = "mean", vjust = -1), 
    #          color = "#D55E00", check_overlap=T )
    ggsave(plot = waterfall_plot, file = paste0(filename_prefix,'_delta_barplot.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300)
  }

waterfall_plotting_numrepeats <- function(distance_data, filename_prefix){
  distance_data$index <- row.names(distance_data)
  distance_mean <- mean(distance_data$delta)
  waterfall_plot<- ggplot(distance_data, aes(x=reorder(index,num_repeats), y=num_repeats)) +
    geom_bar(stat="identity", 
             fill="#0072B2",
             width=0.7,
             position = position_dodge(width=0.4))+
    scale_y_continuous(breaks = round(seq(min(distance_data$distance), max(distance_data$distance), by = 1000),1))+
    labs(title = "Molecule Repeat expansion From hg38 Reference", 
         x = "Molecules", 
         y = "Number of repeat expansions (Units)") +
    theme_classic()+
    theme(axis.text.x = element_blank()) +
    theme(legend.position="none") #+
    # geom_hline(yintercept=75, linetype="dashed", color = "#D55E00")+
    # geom_text(aes(25,75,label = "pathogenic repeat threshold", vjust = -1), 
    #          color = "#D55E00", check_overlap=T )
    ggsave(plot = waterfall_plot, file = paste0(filename_prefix,'_numrepeat_barplot.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300)
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
          stat_function(fun=adjusteddnorm, n=1000, args=list(mm=mod$parameters$mean[i], sd=mod$sd[i], prop=mod$parameters$pro[i]), aes(fill=!!colors[i], color=!!colors[i]), geom="area", alpha = .2)+
          ggtitle(paste0("clustering molecule distances between", "\n labels ", startlabelid, "-", endlabelid))+
          xlab("Distance (bp) between Labels of Interest")+scale_color_discrete(guide=FALSE)+ylab("Number of Molecules")+
          # geom_vline(xintercept=27686, linetype="dashed", color = "#D55E00")+
          # xlim(25000,30000)+
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
    write.csv(simple_result, paste0(out_handle, "_auto_clustering_cluster_simple_result.csv"))
    # completeframe$ClusterNum <- mod$classification
    # cluster_summary <-completeframe[,c("ClusterNum", "distance")]
    # write.csv(mod,
    #                 file=paste0(out_handle, "mod_cluster_summary_noOutlier_delta.csv"),
    #           row.names=FALSE,
    #           quote=FALSE)
    
    return(result)
}

# main function
molecule_distance_delta_ref_main(input_filepath, startlabelid, endlabelid, out_handle)